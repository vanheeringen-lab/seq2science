import os
import re
import time
import pickle
import subprocess
import pandas as pd
from multiprocessing.pool import ThreadPool
from snakemake.utils import validate
from snakemake.logging import logger


# read the samples file
samples = pd.read_csv(config["samples"], sep='\t')
for schema in schemas:
    validate(samples, schema=f"../schemas/samples/{schema}.schema.yaml")
samples['sample'] = samples['sample'].str.strip() 
samples = samples.set_index('sample')
samples.index = samples.index.map(str)


# apply workflow specific changes
if 'assembly' in samples:
    config['assemblies'] = set(samples['assembly'])

if config.get('peak_caller', False):
    config['peak_caller'] = {k: v for d in config['peak_caller'] for k, v in d.items()}

# cut off trailing slashes and make absolute path
for key, value in config.items():
    if '_dir' in key:
        if key in ['result_dir', 'genome_dir', 'rule_dir']:
            value = os.path.abspath(value)
        config[key] = re.split("\/$", value)[0]

# check if paired-end filename suffixes are lexicographically ordered
config['fqext'] = [config['fqext1'], config['fqext2']]
assert sorted(config['fqext'])[0] == config['fqext1']


# check if a sample is single-end or paired end, and store it
logger.info("Checking if samples are single-end or paired-end...")
layout_cachefile = './.snakemake/layouts.p'

def get_layout(sample):
    """ sends a request to ncbi checking whether a sample is single-end or paired-end """
    api_key = config.get('ncbi_key', "")
    if api_key is not "":
        api_key = f'-api_key {api_key}'

    return subprocess.check_output(
        f'''esearch {api_key} -db sra -query {sample} | efetch {api_key} | grep -Po "(?<=<LIBRARY_LAYOUT><)[^/><]*"''',
        shell=True).decode('ascii').rstrip()


# try to load the layout cache, otherwise defaults to empty dictionary
try:
    layout_cache = pickle.load(open(layout_cachefile, "rb"))
except FileNotFoundError:
    layout_cache = {}


tp = ThreadPool(config['ncbi_requests'] // 2)
config['layout'] = {}

# now do a request for each sample that was not in the cache
for sample in [sample for sample in samples.index if sample not in layout_cache]:
    if os.path.exists(expand(f'{{result_dir}}/{{fastq_dir}}/{sample}.{{fqsuffix}}.gz', **config)[0]):
        config['layout'][sample] ='SINGLE'
    elif all(os.path.exists(path) for path in expand(f'{{result_dir}}/{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
        config['layout'][sample] ='PAIRED'
    elif sample.startswith(('GSM', 'SRR', 'ERR', 'DRR')):
        config['layout'][sample] = tp.apply_async(get_layout, (sample,))
        # sleep 1.25 times the minimum required sleep time
        time.sleep(1.25 / (config['ncbi_requests'] // 2))
    else:
        raise ValueError(f"\nsample {sample} was not found..\n"
                         f"We checked for SE file:\n"
                         f"\t{config['result_dir']}/{config['fastq_dir']}/{sample}.{config['fqsuffix']}.gz \n"
                         f"and for PE files:\n"
                         f"\t{config['result_dir']}/{config['fastq_dir']}/{sample}_{config['fqext1']}.{config['fqsuffix']}.gz \n"
                         f"\t{config['result_dir']}/{config['fastq_dir']}/{sample}_{config['fqext2']}.{config['fqsuffix']}.gz \n"
                         f"and since the sample did not start with either GSM, SRR, ERR, and DRR we couldn't find it online..\n")

# now parse the output and store the cache, the local files' layout, and the ones that were fetched online
config['layout'] = {**layout_cache,
                    **{k: (v if isinstance(v, str) else v.get()) for k, v in config['layout'].items()}}

# if new samples were added, update the cache
if len([sample for sample in samples.index if sample not in layout_cache]) is not 0:
    pickle.dump(config['layout'], open(layout_cachefile, "wb"))

logger.info("Done!\n\n")


# after all is done, log (print) the configuration
logger.info("CONFIGURATION VARIABLES:")
for key, value in config.items():
     logger.info(f"{key: <23}: {value}")
logger.info("\n\n")


# If hmmratac peak caller, check if all samples are paired-end
if config.get('peak_caller', False) and 'hmmratac' in config['peak_caller']:
    assert all([config['layout'][sample] == 'PAIRED' for sample in samples.index]), \
    "HMMRATAC requires all samples to be paired end"


# Find conda directories. Does not work with singularity.
def conda_path(yaml):
    """ Find the path to a conda directory """
    import hashlib
    import os.path

    env_file = os.path.abspath(yaml)
    env_dir = os.path.join(os.getcwd(), ".snakemake", "conda")

    md5hash = hashlib.md5()
    md5hash.update(env_dir.encode())
    with open(env_file, 'rb') as f:
        content = f.read()
    md5hash.update(content)
    hash = md5hash.hexdigest()[:8]
    path = os.path.join(env_dir, hash)
    return path
