import os.path
import re
import time
import math
import norns
import psutil
import pickle
import subprocess
import numpy as np
import pandas as pd
from multiprocessing.pool import ThreadPool
from snakemake.utils import validate
from snakemake.logging import logger


# make sure that difficult data-types (yaml objects) are in correct data format
for kw in ['aligner', 'bam_sorter']:
    if isinstance(config.get(kw, None), str):
        config[kw] = {config[kw]: {}}

# validate and complement the config dict
for schema in config_schemas:
    validate(config, schema=f"../schemas/config/{schema}.schema.yaml")

# check if paired-end filename suffixes are lexicographically ordered
config['fqext'] = [config['fqext1'], config['fqext2']]
assert sorted(config['fqext'])[0] == config['fqext1']

# apply workflow specific changes
# for atac-seq
if config.get('peak_caller', False):
    config['peak_caller'] = {k: v for k,v in config['peak_caller'].items()}

    # if hmmratac peak caller, check if all samples are paired-end
    if 'hmmratac' in config['peak_caller']:
        assert all([config['layout'][sample] == 'PAIRED' for sample in samples.index]), \
        "HMMRATAC requires all samples to be paired end"

# for rna-seq
for conf_dict in ['aligner', 'diffexp']:
    if config.get(conf_dict, False):
        dict_key = list(config[conf_dict].keys())[0]
        for k, v in list(config[conf_dict].values())[0].items():
            config[k] = v
        config[conf_dict] = dict_key

# for alignment
if config.get('bam_sorter', False):
    config['bam_sort_order'] = list(config['bam_sorter'].values())[0]
    config['bam_sorter'] = list(config['bam_sorter'].keys())[0]


# read and validate the samples file
samples = pd.read_csv(config["samples"], sep='\t')
for schema in sample_schemas:
    validate(samples, schema=f"../schemas/samples/{schema}.schema.yaml")
samples['sample'] = samples['sample'].str.strip()
samples = samples.set_index('sample')
samples.index = samples.index.map(str)

# make sure that our samples.tsv and configuration work together (e.g. conditions)
if 'condition' in samples:
    if 'hmmratac' in config['peak_caller']:
        assert config['combine_replicates'] == 'idr', \
        f'HMMRATAC peaks can only be combined through idr'

    for condition in set(samples['condition']):
        for assembly in set(samples[samples['condition'] == condition]['assembly']):
            nr_samples = len(samples[(samples['condition'] == condition) & (samples['assembly'] == assembly)])
            assert nr_samples >= 2,\
            f'When specifying conditions every condition needs at least two samples, however you gave {nr_samples}'\
            f' sample for condition {condition} and assembly {assembly}'

            if config.get('combine_replicates', '') == 'idr':
                assert nr_samples == 2,\
                f'For IDR to work you need two samples per condition, however you gave {nr_samples} samples for'\
                f' condition {condition} and assembly {assembly}'


try:
    user_config = norns.config(config_file=workflow.overwrite_configfile)
except:
    user_config = norns.config(config_file='config.yaml')

# make absolute paths, cut off trailing slashes
for key, value in config.items():
    if '_dir' in key:
        if key in ['result_dir', 'genome_dir', 'rule_dir'] or key in user_config:
            value = os.path.abspath(value)
        config[key] = re.split("\/$", value)[0]


# check if a sample is single-end or paired end, and store it
logger.info("Checking if samples are single-end or paired-end...")
layout_cachefile = './.snakemake/layouts.p'

def get_layout(sample):
    """ sends a request to ncbi checking whether a sample is single-end or paired-end """
    api_key = config.get('ncbi_key', "")
    if api_key is not "":
        api_key = f'-api_key {api_key}'

    try:
        layout = subprocess.check_output(
            f'''esearch {api_key} -db sra -query {sample} | efetch {api_key} | grep -Po "(?<=<LIBRARY_LAYOUT><)[^ /><]*"''',
            shell=True).decode('ascii').rstrip()
        if layout not in ['PAIRED', 'SINGLE']:
            raise ValueError(f"Sample {sample} was found to be {layout}-end, however the only"
                             f"acceptable library layouts are SINGLE-end and PAIRED-end.")
        return layout
    except subprocess.CalledProcessError:
        raise ValueError(f"The command to lookup sample {sample} online failed!\n"
                         f"Are you sure this sample exists..? Downloading samples with restricted "
                         f"access is currently not supported. We advise you to download the sample "
                         f"manually, and continue the pipeline from there on.")


# try to load the layout cache, otherwise defaults to empty dictionary
try:
    layout_cache = pickle.load(open(layout_cachefile, "rb"))
except FileNotFoundError:
    layout_cache = {}


tp = ThreadPool(config.get('ncbi_requests', 3) // 2)
config['layout'] = {}

# now do a request for each sample that was not in the cache
for sample in [sample for sample in samples.index if sample not in layout_cache]:
    if os.path.exists(expand(f'{{fastq_dir}}/{sample}.{{fqsuffix}}.gz', **config)[0]):
        config['layout'][sample] ='SINGLE'
    elif all(os.path.exists(path) for path in expand(f'{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
        config['layout'][sample] ='PAIRED'
    elif sample.startswith(('GSM', 'SRR', 'ERR', 'DRR')):
        config['layout'][sample] = tp.apply_async(get_layout, (sample,))
        # sleep 1.25 times the minimum required sleep time
        time.sleep(1.25 / (config.get('ncbi_requests', 3) // 2))
    else:
        raise ValueError(f"\nsample {sample} was not found..\n"
                         f"We checked for SE file:\n"
                         f"\t/{config['fastq_dir']}/{sample}.{config['fqsuffix']}.gz \n"
                         f"and for PE files:\n"
                         f"\t{config['fastq_dir']}/{sample}_{config['fqext1']}.{config['fqsuffix']}.gz \n"
                         f"\t{config['fastq_dir']}/{sample}_{config['fqext2']}.{config['fqsuffix']}.gz \n"
                         f"and since the sample did not start with either GSM, SRR, ERR, and DRR we couldn't find it online..\n")

# now parse the output and store the cache, the local files' layout, and the ones that were fetched online
config['layout'] = {**layout_cache,
                    **{k: (v if isinstance(v, str) else v.get()) for k, v in config['layout'].items()}}

assert all(layout in ['SINGLE', 'PAIRED'] for sample, layout in config['layout'].items())

# if new samples were added, update the cache
if len([sample for sample in samples.index if sample not in layout_cache]) is not 0:
    pickle.dump(config['layout'], open(layout_cachefile, "wb"))

logger.info("Done!\n\n")

# if samples are merged add the layout of the condition to the config
if 'condition' in samples and config.get('combine_replicates', "") == 'merge':
    for sample in samples.index:
        condition = samples.loc[sample, 'condition']
        config['layout'][condition] = config['layout'][sample]

# after all is done, log (print) the configuration
logger.info("CONFIGURATION VARIABLES:")
for key, value in config.items():
     logger.info(f"{key: <23}: {value}")
logger.info("\n\n")


# if differential gene expression analysis is used, check all contrasts
if config.get('contrasts', False):
    old_contrasts = list(config["contrasts"])
    for contrast in old_contrasts:
        original_contrast = contrast

        # remove whitespace
        contrast = contrast.replace(" ", "")
        contrast = contrast.replace("~", "")

        # remove batch
        batch = None
        if '+' in contrast:
            batch = contrast.split('+')[0]
            contrast = contrast.split('+')[1]

        # parse contrast
        contrast = contrast.split('_')

        # Check if the contrast can be recognized
        assert contrast[0] in samples.columns, \
            f'In contrast design {original_contrast}, {contrast[0]} does not match any column name in {config["samples"]}'
        if batch is not None:
            assert batch in samples.columns, \
                f'In contrast design {original_contrast}, the batch effect {batch} does not match any column name in {config["samples"]}'

        l = len(contrast)
        if l == 1:
            # check if contrast column has exactly 2 factor levels (per assembly)
            tmp = samples[['assembly', contrast[0]]].dropna()
            factors = pd.DataFrame(tmp.groupby('assembly')[contrast[0]].nunique())
            assert all(factors[contrast[0]] == 2),\
                f'Your contrast design, {original_contrast}, contains only a column name ({contrast[0]}), '\
                f'If you wish to compare all groups in this column, add a reference group. \n'\
                f'number of groups found (per assembly): \n'\
                f'{factors[contrast[0]]}'
        if l > 1:
            # check if contrast column contains the groups
            for group in contrast[1:]:
                if group != 'all':
                    assert str(group) in str(samples[contrast[0]].tolist()),\
                    f'Your contrast design contains group {group}, '
                    f'which cannot be found in column {contrast[0]} of {config["samples"]}'


# regex compatible string of all elements in the samples.tsv column given by the input
def any_given(*args):
    elements = []
    for column_name in args:
        if column_name in samples:
            elements.extend(samples[column_name].unique())
        elif column_name is 'sample':
            elements.extend(samples.index.unique())

    return '|'.join(elements)

# set global constraints on wildcards ({{sample}} or {{assembly}})
if 'assembly' in samples:
    wildcard_constraints:
        sample=any_given('sample', 'condition'),
        assembly=any_given('assembly')
else:
    wildcard_constraints:
        sample=any_given('sample', 'condition')


# set default parameters (parallel downloads and memory)
def convert_size(size_bytes, order=None):
    # https://stackoverflow.com/questions/5194057/better-way-to-convert-file-sizes-in-python/14822210#14822210
    size_name = ["b", "kb", "mb", "gb"]
    if order is None:
        order = int(math.floor(math.log(size_bytes, 1024)))
    s = int(size_bytes // math.pow(1024, order))
    return s, size_name[order]


def add_default_resources(func):
    # https://stackoverflow.com/questions/6200270/decorator-to-print-function-call-details-parameters-names-and-effective-values/6278457#6278457
    def wrapper(*args, **kwargs):
        # by default only one download in parallel (workflow fails on multiple on a single node)
        kwargs['resources'] = {**{'parallel_downloads': 1}, **kwargs['resources']}

        # when the user specifies memory, use this and give a warning if it surpasses local memory
        # (surpassing does not always have to be an issue -> cluster execution)
        # if none specified set the memory to the max available on the computer
        mem = psutil.virtual_memory()
        if kwargs['resources'].get('mem_mb'):
            if kwargs['resources']['mem_mb'] > convert_size(mem.total, 2)[0]:
                logger.info(f"WARNING: The specified ram ({kwargs['resources']['mem_mb']} mb) surpasses the local machine\'s RAM ({convert_size(mem.total, 2)[0]} mb)")

            kwargs['resources'] = {**kwargs['resources'],
                                   **{'mem_mb': np.clip(kwargs['resources']['mem_mb'], 0, convert_size(mem.total, 2)[0])}}
        else:
            if kwargs['resources'].get('mem_gb', 0) > convert_size(mem.total, 3)[0]:
                logger.info(f"WARNING: The specified ram ({kwargs['resources']['mem_gb']} gb) surpasses the local machine\'s RAM ({convert_size(mem.total, 3)[0]} gb)")

            kwargs['resources'] = {**kwargs['resources'],
                                   **{'mem_gb': np.clip(kwargs['resources'].get('mem_gb', 9999), 0, convert_size(mem.total, 3)[0])}}

        return func(*args, **kwargs)
    return wrapper


# now add the wrapper to the workflow execute function
workflow.execute = add_default_resources(workflow.execute)



# functional but currently unused
# # find conda directories. Does not work with singularity.
# def conda_path(yaml):
#     """ Find the path to a conda directory """
#     import hashlib
#     import os.path
#
#     env_file = os.path.abspath(yaml)
#     env_dir = os.path.join(os.getcwd(), ".snakemake", "conda")
#
#     md5hash = hashlib.md5()
#     md5hash.update(env_dir.encode())
#     with open(env_file, 'rb') as f:
#         content = f.read()
#     md5hash.update(content)
#     dir_hash = md5hash.hexdigest()[:8]
#     path = os.path.join(env_dir, dir_hash)
#     return path
