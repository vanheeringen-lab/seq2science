import os.path
import glob
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


# check config file for correct directory names
for key, value in config.items():
    if '_dir' in key:
        # allow tilde as the first character only
        assert not re.match('[^A-Za-z0-9_.\-/~]+', value[0]) and not re.match('[^A-Za-z0-9_.\-/]+', value[1:]) and not " " in value, \
            ("\n" + "In the config.yaml you set '" + key + "' to '" + value + "'. Please use file paths that only contain letters, " +
            "numbers and any of the following symbols: underscores (_), periods (.), and minuses (-) (the first symbol may be a tilde (~)).\n")
        config[key] = os.path.expanduser(value)

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

# read and sanitize the samples file
samples = pd.read_csv(config["samples"], sep='\t')
# sanitize column names
samples.columns = samples.columns.str.strip()
assert all([col[0:7] not in ["Unnamed", ''] for col in samples]), \
    ("\nEncountered unnamed column in " + config["samples"] +
     ".\nColumn names: " + str(', '.join(samples.columns)) + '.\n')
assert not any(samples.columns.str.contains('[^A-Za-z0-9_.\-]+', regex=True)), \
    ("\n" + config["samples"] + " may only contain letters, numbers and " +
    "underscores (_), periods (.), or minuses (-).\n")
# sanitize table content
samples = samples.applymap(lambda x: str(x).strip())
assert not any([any(samples[col].str.contains('[^A-Za-z0-9_.\-]+', regex=True)) for col in samples]), \
    ("\n" + config["samples"] + " may only contain letters, numbers and " +
    "underscores (_), periods (.), or minuses (-).\n")
assert len(samples["sample"]) == len(set(samples["sample"])), \
    ("\nDuplicate samples found in " + config["samples"] + ":\n" +
     samples[samples.duplicated(['sample'], keep=False)].to_string() + '\n')

# validate samples file
for schema in sample_schemas:
    validate(samples, schema=f"../schemas/samples/{schema}.schema.yaml")
samples = samples.set_index('sample')
samples.index = samples.index.map(str)


# apply workflow specific changes...
# ...for atac-seq
if config.get('peak_caller', False):
    config['peak_caller'] = {k: v for k,v in config['peak_caller'].items()}

    # if hmmratac peak caller, check if all samples are paired-end
    if 'hmmratac' in config['peak_caller']:
        assert all([config['layout'][sample] == 'PAIRED' for sample in samples.index]), \
        "HMMRATAC requires all samples to be paired end"

# ...for alignment and rna-seq
for conf_dict in ['aligner', 'diffexp']:
    if config.get(conf_dict, False):
        dict_key = list(config[conf_dict].keys())[0]
        for k, v in list(config[conf_dict].values())[0].items():
            config[k] = v
        config[conf_dict] = dict_key

# ...for alignment
if config.get('bam_sorter', False):
    config['bam_sort_order'] = list(config['bam_sorter'].values())[0]
    config['bam_sorter'] = list(config['bam_sorter'].keys())[0]

try:
    user_config = norns.config(config_file=workflow.overwrite_configfile)
except:
    user_config = norns.config(config_file='config.yaml')


# make sure that our samples.tsv and configuration work together...
# ...on replicates
if 'condition' in samples:
    if config['combine_replicates'] != 'merge':
        if 'hmmratac' in config['peak_caller']:
            assert config['combine_replicates'] == 'idr', \
            f'HMMRATAC peaks can only be combined through idr'

    for condition in set(samples['condition']):
        for assembly in set(samples[samples['condition'] == condition]['assembly']):
            nr_samples = len(samples[(samples['condition'] == condition) & (samples['assembly'] == assembly)])

            if config.get('combine_replicates', '') == 'idr':
                assert nr_samples == 2,\
                f'For IDR to work you need two samples per condition, however you gave {nr_samples} samples for'\
                f' condition {condition} and assembly {assembly}'

# ...on DE contrasts
def parse_DE_contrasts(de_contrast):
    """
    Extract batch and contrast groups from a DE contrast design
    """
    original_contrast = de_contrast

    # remove whitespaces (and '~'s if used)
    de_contrast = de_contrast.replace(" ", "").replace("~", "")

    # split and store batch effect
    batch = None
    if '+' in de_contrast:
        batch =  de_contrast.split('+')[0]
        de_contrast = de_contrast.split('+')[1]

    # parse contrast
    parsed_contrast = de_contrast.split('_')
    return original_contrast, parsed_contrast, batch

if config.get('contrasts', False):
    # check differential gene expression contrasts
    old_contrasts = list(config["contrasts"])
    for contrast in old_contrasts:
        original_contrast, parsed_contrast, batch = parse_DE_contrasts(contrast)

        # Check if the column names can be recognized in the contrast
        assert parsed_contrast[0] in samples.columns and parsed_contrast[0] not in ["sample", "assembly"], \
            ('\nIn contrast design "' + original_contrast + '", "' + parsed_contrast[0] +
             '" does not match any valid column name in ' + config["samples"] + '.\n')
        if batch is not None:
            assert batch in samples.columns and batch not in ["sample", "assembly"], \
                ('\nIn contrast design "' + original_contrast + '", the batch effect "' +
                 batch + '" does not match any valid column name in ' + config["samples"] + '.\n')

        # Check if the groups described by the contrast can be identified and found in samples.tsv
        l = len(parsed_contrast)
        assert l < 4, ("\nA differential expression contrast couldn't be parsed correctly. "
                       + str(l-1) + " groups were found in \n" + original_contrast + ' (groups: ' +
                       ', '.join(parsed_contrast[1:]) + ').' +
                       '\nPlease do not use whitespaces or underscores in your contrast, ' +
                       '\nor in the columns in ' + config['samples'] + ' referenced by your contrast.\n')
        if l == 1:
            # check if contrast column has exactly 2 factor levels (per assembly)
            tmp = samples[['assembly', parsed_contrast[0]]].dropna()
            factors = pd.DataFrame(tmp.groupby('assembly')[parsed_contrast[0]].nunique())
            assert all(factors[parsed_contrast[0]] == 2),\
                ('\nYour contrast design, ' + original_contrast +
                 ', contains only a column name (' + parsed_contrast[0] +
                 '). \nIf you wish to compare all groups in this column, add a reference group. ' +
                 'Number of groups found (per assembly): \n' + str(factors[parsed_contrast[0]]))
        else:
            # check if contrast column contains the groups
            for group in parsed_contrast[1:]:
                if group != 'all':
                    assert str(group) in [str(i) for i in samples[parsed_contrast[0]].tolist()],\
                        ('\nYour contrast design contains group ' + group +
                        ' which cannot be found in column ' + parsed_contrast[0] +
                         'of' + config["samples"] + '.\n')


# make absolute paths, cut off trailing slashes
for key, value in config.items():
    if '_dir' in key:
        if key in ['result_dir', 'genome_dir', 'rule_dir'] or key in user_config:
            value = os.path.abspath(value)
        config[key] = re.split("\/$", value)[0]

# nest default directories in result_dir
for key, value in config.items():
    if '_dir' in key:
        if key not in ['result_dir', 'genome_dir', 'rule_dir'] and key not in user_config:
            config[key] = os.path.join(config['result_dir'], config[key])


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

# now only keep the layout of samples that are in samples.tsv
config['layout'] = {key: value for key, value in config['layout'].items() if key in samples.index}

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


def any_given(*args):
    """
    returns a regex compatible string of all elements in the samples.tsv column given by the input
    """
    elements = []
    for column_name in args:
        if column_name in samples:
            elements.extend(samples[column_name])
        elif column_name is 'sample':
            elements.extend(samples.index)

    return '|'.join(set(elements))

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
# workflow.execute = add_default_resources(workflow.execute)


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
