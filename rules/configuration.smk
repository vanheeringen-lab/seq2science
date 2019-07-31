import os
import re
import time
import genomepy
import subprocess
import pandas as pd
from multiprocessing.pool import ThreadPool
from snakemake.utils import validate
from snakemake.logging import logger


# read the samples file
samples = pd.read_csv(config["samples"], sep='\t')
validate(samples, schema=f"../schemas/{schema}")
samples = samples.set_index('sample')
samples.index = samples.index.map(str)


# TODO: maybe move to other locations? Cleaner?
# apply workflow specific changes
if 'condition' in samples:
    samples['condition'] = samples['condition'].str.replace(" ","")

if 'assembly' in samples:
    config['assemblies'] = set(samples['assembly'])

if config['peak_caller']:
    config['peak_caller'] = {k: v for d in config['peak_caller'] for k, v in d.items()}
    assert all(key in ['macs2', 'genrich'] for key in config['peak_caller'].keys())

# cut off trailing slashes
for path in ['result_dir', 'genome_dir', 'log_dir']:
    config[path] = re.split("\/$", config[path])[0]

# check if paired-end filename suffixes are lexicograpically ordered
config['fqext'] = [config['fqext1'], config['fqext2']]
assert sorted(config['fqext'])[0] == config['fqext1']


# check if a sample is single-end or paired end, and store it
def get_layout(sample):
    """ sends a request to ncbi checking whether a sample is single-end or paired-end """
    api_key = config.get('ncbi_key', "")
    if api_key is not "":
        api_key = f'-api_key {api_key}'

    return sample, subprocess.check_output(
        f'''esearch {api_key} -db sra -query {sample} | efetch {api_key} | grep -Po "(?<=<LIBRARY_LAYOUT><)[^/><]*"''',
        shell=True).decode('ascii').rstrip()


results = []
tp = ThreadPool(config['ncbi_requests'] // 2)
config['layout'] = {}

# now do a request for each sample
for sample in samples.index:
    if   os.path.exists(expand(f'{{result_dir}}/{{fastq_dir}}/SE/{sample}.{{fqsuffix}}.gz', **config)[0]):
        config['layout'][sample] ='SINGLE'
    elif all(os.path.exists(path) for path in expand(f'{{result_dir}}/{{fastq_dir}}/PE/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
        config['layout'][sample] ='PAIRED'
    else:
        results.append(tp.apply_async(get_layout, (sample,)))
        # sleep 1.25 times the minimum required sleep time
        time.sleep(1.25 / (config['ncbi_requests'] // 2))

# now parse the output and store in config
config['layout'] = {**config['layout'], **{r.get()[0]: r.get()[1] for r in results}}


# Do onstart/onexit things
onstart:
    # turn off genomepy plugins, since they slow down the 'downloading' process
    if 'genomepy' in sys.modules:
        # get the genomepy settings
        config['active_plugins'] = [p.name() for p in genomepy.plugin.get_active_plugins()]

        # disable genomepy plugins
        for plugin in ['bowtie2', 'bwa', 'gmap', 'hisat2', 'minimap2']:
            genomepy.plugin.deactivate(plugin)


def onexit(config):
    if 'genomepy' in sys.modules:
        for plugin in genomepy.plugin.plugins:
            if plugin in config['active_plugins']:
                genomepy.plugin.activate(plugin)


onerror:
    onexit(config)


onsuccess:
    onexit(config)


# after all is done, log (print) the configuration
logger.info("CONFIGURATION VARIABLES:")
for key, value in config.items():
     logger.info(f"{key: <23}: {value}")
logger.info("\n\n")
