import os
import re
import time
import subprocess
import pandas as pd
from ast import literal_eval
from hashlib import blake2b
from multiprocessing.pool import ThreadPool
from snakemake.utils import validate
from snakemake.logging import logger


# read the samples file
samples = pd.read_csv(config["samples"], sep='\t')
validate(samples, schema=f"../schemas/{schema}")
samples['sample'] = samples['sample'].str.strip()
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

# cut off trailing slashes and make absolute path
for key, value in config.items():
    if '_dir' in key:
        if key in ['result_dir', 'genome_dir']:
            value = os.path.abspath(value)
        config[key] = re.split("\/$", value)[0]

# check if paired-end filename suffixes are lexicograpically ordered
config['fqext'] = [config['fqext1'], config['fqext2']]
assert sorted(config['fqext'])[0] == config['fqext1']


# check if a sample is single-end or paired end, and store it
logger.info("Checking if samples are single-end or paired-end...")
# layout file with a unique name depending on the samples contained
identifier = str(blake2b(''.join(samples.sort_index().index).encode(), digest_size=12).hexdigest())
layout_file = './.snakemake/layouts_' + identifier
if not os.path.exists(layout_file):
    def get_layout(sample):
        """ sends a request to ncbi checking whether a sample is single-end or paired-end """
        api_key = config.get('ncbi_key', "")
        if api_key is not "":
            api_key = f'-api_key {api_key}'

        return sample, subprocess.check_output(
            f'''esearch {api_key} -db sra -query {sample} | efetch {api_key} | grep -Po "(?<=<LIBRARY_LAYOUT><)[^/><]*"''',
            shell=True).decode('ascii').rstrip()


    tp = ThreadPool(config['ncbi_requests'] // 2)
    layout = {}

    # now do a request for each sample
    for sample in samples.index:
        if os.path.exists(expand(f'{{result_dir}}/{{fastq_dir}}/{sample}.{{fqsuffix}}.gz', **config)[0]):
            layout[sample] = 'SINGLE'
        elif all(os.path.exists(path) for path in expand(f'{{result_dir}}/{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
            layout[sample] = 'PAIRED'
        else:
            layout[sample] = get_layout(sample)[1]
            # sleep 1.25 times the minimum required sleep time
            time.sleep(1.25 / (config['ncbi_requests'] // 2))

    # save layout to file
    with open(layout_file, 'w') as f:
        f.write(str(layout))
else:
    with open(layout_file, 'r') as f:
        layout = literal_eval(f.read())

# now store in config
config['layout'] = layout
logger.info("Done!\n\n")


# after all is done, log (print) the configuration
logger.info("CONFIGURATION VARIABLES:")
for key, value in config.items():
     logger.info(f"{key: <23}: {value}")
logger.info("\n\n")