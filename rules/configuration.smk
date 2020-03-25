import math
import os.path
import psutil
import pickle
import re
import shutil
import subprocess
import time

import norns
import numpy as np
import pandas as pd
import urllib.request
from bs4 import BeautifulSoup
from multiprocessing.pool import ThreadPool
from filelock import FileLock

from snakemake.logging import logger
from snakemake.utils import validate
from snakemake.utils import min_version


# make sure the snakemake version corresponds to version in environment
min_version("5.10")

# check config file for correct directory names
for key, value in config.items():
    if '_dir' in key:
        assert value not in [None, '', ' '], f"\n{key} cannot be empty. For current directory, set '{key}: .'\n"
        # allow tilde as the first character, \w = all letters and numbers
        assert re.match('^[~\w_./-]*$', value[0]) and re.match('^[\w_./-]*$', value[1:]), \
            (f"\nIn the config.yaml you set '{key}' to '{value}'. Please use file paths that only contain letters, " +
            "numbers and any of the following symbols: underscores (_), periods (.), and minuses (-). The first character may also be a tilde (~).\n")
        config[key] = os.path.expanduser(value)

# make sure that difficult data-types (yaml objects) are in correct data format
for kw in ['aligner', 'quantifier', 'bam_sorter']:
    if isinstance(config.get(kw, None), str):
        config[kw] = {config[kw]: {}}

# validate and complement the config dict
for schema in config_schemas:
    validate(config, schema=f"../schemas/config/{schema}.schema.yaml")

# check if paired-end filename suffixes are lexicographically ordered
config['fqext'] = [config['fqext1'], config['fqext2']]
assert sorted(config['fqext'])[0] == config['fqext1'], \
    ("\nThe paired-end filename suffixes must be lexicographically ordered!\n" +
     f"Example suffixes: fqext1: R1, fqext2: R2\n" +
     f"Your suffixes:    fqext1: {config['fqext1']}, fqext2: {config['fqext2']}\n")

# read and sanitize the samples file
samples = pd.read_csv(config["samples"], sep='\t', dtype='str')

# sanitize column names
samples.columns = samples.columns.str.strip()
assert all([col[0:7] not in ["Unnamed", ''] for col in samples]), \
    (f"\nEncountered unnamed column in {config['samples']}.\n" +
     f"Column names: {str(', '.join(samples.columns))}.\n")
assert not any(samples.columns.str.contains('[^A-Za-z0-9_.\-%]+', regex=True)), \
    (f"\n{config['samples']} may only contain letters, numbers and " +
    "percentage signs (%), underscores (_), periods (.), or minuses (-).\n")

# sanitize table content
if 'replicate' in samples:
    samples['replicate'] = samples['replicate'].mask(pd.isnull, samples['sample'])
    # if there is nothing to merge, drop the column. keep it simple
    if samples['replicate'].tolist() == samples['sample'].tolist() or config.get('technical_replicates') == 'keep':
        samples = samples.drop(columns=['replicate'])
if 'condition' in samples:
    samples['condition'] = samples['condition'].mask(pd.isnull, samples['replicate']) if 'replicate' in samples else \
        samples['condition'].mask(pd.isnull, samples['sample'])
    # if there is nothing to merge, drop the column. keep it simple
    if samples['condition'].tolist() == samples['sample'].tolist() or config.get('biological_replicates') == 'keep':
        samples = samples.drop(columns=['condition'])

if 'replicate' in samples:
    r = samples[['assembly', 'replicate']].drop_duplicates().set_index('replicate')
    for replicate in r.index:
        assert len(r[r.index == replicate]) == 1, \
            ("\nReplicate names must be different between assemblies.\n" +
             f"Replicate name '{replicate}' was found in assemblies {r[r.index == replicate]['assembly'].tolist()}.")

assert not any([any(samples[col].str.contains('[^A-Za-z0-9_.\-%]+', regex=True)) for col in samples]), \
    (f"\n{config['samples']} may only contain letters, numbers and " +
    "percentage signs (%), underscores (_), periods (.), or minuses (-).\n")
assert len(samples["sample"]) == len(set(samples["sample"])), \
    (f"\nDuplicate samples found in {config['samples']}:\n" +
     f"{samples[samples.duplicated(['sample'], keep=False)].to_string()}\n")

for idx in samples.index:
    if "condition" in samples:
        assert idx not in samples["condition"].values, f"sample names, conditions, and replicates can not overlap. " \
                                                       f"Sample {idx} can not also occur as a condition"
    if "replicate" in samples:
        assert idx not in samples["replicate"].values, f"sample names, conditions, and replicates can not overlap. " \
                                                       f"Sample {idx} can not also occur as a replicate"

if "condition" in samples and "replicate" in samples:
    for cond in samples["condition"]:
        assert cond not in samples["replicate"].values, f"sample names, conditions, and replicates can not overlap. " \
                                                        f"Condition {cond} can not also occur as a replicate"

# validate samples file
for schema in sample_schemas:
    validate(samples, schema=f"../schemas/samples/{schema}.schema.yaml")
samples = samples.set_index('sample')
samples.index = samples.index.map(str)


# apply workflow specific changes...
# ...for atac-seq
if config.get('peak_caller', False):
    config['peak_caller'] = {k: v for k,v in config['peak_caller'].items()}

    # if genrich is peak caller, make sure to not double shift reads
    if 'genrich' in config['peak_caller']:
        # always turn of genrich shift, since we handle that with deeptools
        if '-j' in config['peak_caller']['genrich'] and not '-D' in config['peak_caller']['genrich']:
            config['peak_caller']['genrich'] += ' -D'

    # if hmmratac peak caller, check if all samples are paired-end
    if 'hmmratac' in config['peak_caller']:
        assert all([config['layout'][sample] == 'PAIRED' for sample in samples.index]), \
        "HMMRATAC requires all samples to be paired end"

    config['macs2_types'] = ['control_lambda.bdg', 'peaks.xls', 'treat_pileup.bdg']
    if 'macs2' in config['peak_caller'] and '--broad' in config['peak_caller']['macs2']:
        config['macs2_types'].extend(['peaks.broadPeak', 'peaks.gappedPeak'])
    else:
        config['macs2_types'].extend(['summits.bed', 'peaks.narrowPeak'])


# ...for alignment and rna-seq
for conf_dict in ['aligner', 'quantifier', 'diffexp']:
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
    user_config = norns.config(config_file=workflow.overwrite_configfiles[-1])
except:
    user_config = norns.config(config_file='config.yaml')


# make sure that our samples.tsv and configuration work together...
# ...on biological replicates
if 'condition' in samples and config.get('biological_replicates', '') != 'keep':
    if 'hmmratac' in config['peak_caller']:
        assert config.get('biological_replicates', '') == 'idr', \
        f'HMMRATAC peaks can only be combined through idr'

    for condition in set(samples['condition']):
        for assembly in set(samples[samples['condition'] == condition]['assembly']):
            if 'replicate' in samples and config.get('technical_replicates') == 'merge':
                nr_samples = len(set(samples[(samples['condition'] == condition) & (samples['assembly'] == assembly)]['replicate']))
            else:
                nr_samples = len(samples[(samples['condition'] == condition) & (samples['assembly'] == assembly)])

            if config.get('biological_replicates', '') == 'idr':
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
            (f'\nIn contrast design "{original_contrast}", "{parsed_contrast[0]} ' +
             f'does not match any valid column name in {config["samples"]}.\n')
        if batch is not None:
            assert batch in samples.columns and batch not in ["sample", "assembly"], \
                (f'\nIn contrast design "{original_contrast}", the batch effect "{batch}" ' +
                 f'does not match any valid column name in {config["samples"]}.\n')

        # Check if the groups described by the contrast can be identified and found in samples.tsv
        l = len(parsed_contrast)
        assert l < 4, ("\nA differential expression contrast couldn't be parsed correctly.\n" +
                       f"{str(l-1)} groups were found in '{original_contrast}' " +
                       f"(groups: {', '.join(parsed_contrast[1:])}).\n\n" +
                       f'Tip: do not use underscores in the columns of {config["samples"]} referenced by your contrast.\n')
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
                         ' of ' + config["samples"] + '.\n')


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
layout_cachefile = os.path.expanduser('~/.config/snakemake/layouts.p')
layout_cachefile_lock = os.path.expanduser('~/.config/snakemake/layouts.p.lock')

def get_layout_eutils(sample):
    """
    Sends a request to ncbi checking whether a sample is single-end or paired-end.
    Robust method (always returns result), however manually filled in by uploader, so can be wrong.
    Complementary method of get_layout_trace
    """
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


def get_layout_trace(sample):
    """
    Parse the ncbi trace website to check if a read has 1, 2, or 3 spots.
    Will fail if sample is not on ncbi database, however does not have the problem that uploader
    filled out the form wrongly.
    Complementary method of get_layout_eutils
    """
    try:
        url = f"https://www.ncbi.nlm.nih.gov/sra/?term={sample}"

        conn = urllib.request.urlopen(url)
        html = conn.read()

        soup = BeautifulSoup(html, features="html5lib")
        links = soup.find_all('a')

        for tag in links:
            link = tag.get('href', None)
            if link is not None and 'SRR' in link:
                trace_conn = urllib.request.urlopen("https:" + link)
                trace_html = trace_conn.read()
                x = re.search("This run has (\d) read", str(trace_html))

                # if there are spots without info, then just ignore this sample
                if len(re.findall(", average length: 0", str(trace_html))) > 0:
                    break
                elif x.group(1) == '1':
                    return 'SINGLE'
                elif x.group(1) == '2':
                    return 'PAIRED'
    except:
        pass
    return None

# do this locked to avoid parallel ncbi requests with the same key, and to avoid
# multiple writes/reads at the same time to layouts.p
if not os.path.exists(os.path.dirname(layout_cachefile_lock)):
    os.makedirs(os.path.dirname(layout_cachefile_lock))

# let's ignore locks that are older than 5 minutes
if os.path.exists(layout_cachefile_lock) and \
        time.time() - os.stat(layout_cachefile_lock).st_mtime > 5 * 60:
    os.remove(layout_cachefile_lock)

with FileLock(layout_cachefile_lock):
    # try to load the layout cache, otherwise defaults to empty dictionary
    try:
        layout_cache = pickle.load(open(layout_cachefile, "rb"))
    except FileNotFoundError:
        layout_cache = {}


    trace_tp = ThreadPool(20)
    eutils_tp = ThreadPool(config.get('ncbi_requests', 3) // 2)

    trace_layout = {}
    config['layout'] = {}

    # now do a request for each sample that was not in the cache
    for sample in [sample for sample in samples.index if sample not in layout_cache]:
        config['layout'][sample] = eutils_tp.apply_async(get_layout_eutils, (sample,))
        if os.path.exists(expand(f'{{fastq_dir}}/{sample}.{{fqsuffix}}.gz', **config)[0]):
            config['layout'][sample] ='SINGLE'
        elif all(os.path.exists(path) for path in expand(f'{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
            config['layout'][sample] ='PAIRED'
        elif sample.startswith(('GSM', 'SRR', 'ERR', 'DRR')):
            # config['layout'][sample] = eutils_tp.apply_async(get_layout_eutils, (sample,))
            trace_layout[sample] = trace_tp.apply_async(get_layout_trace, (sample,))

            # sleep 1.25 times the minimum required sleep time so eutils don't complain
            time.sleep(1.25 / (config.get('ncbi_requests', 3) // 2))
        else:
            raise ValueError(f"\nsample {sample} was not found..\n"
                             f"We checked for SE file:\n"
                             f"\t{config['fastq_dir']}/{sample}.{config['fqsuffix']}.gz \n"
                             f"and for PE files:\n"
                             f"\t{config['fastq_dir']}/{sample}_{config['fqext1']}.{config['fqsuffix']}.gz \n"
                             f"\t{config['fastq_dir']}/{sample}_{config['fqext2']}.{config['fqsuffix']}.gz \n"
                             f"and since the sample did not start with either GSM, SRR, ERR, and DRR we couldn't find it online..\n")

    # now parse the output and store the cache, the local files' layout, and the ones that were fetched online
    config['layout'] = {**layout_cache,
                        **{k: (v if isinstance(v, str) else v.get()) for k, v in config['layout'].items()},
                        **{k: v.get() for k, v in trace_layout.items() if v.get() is not None}}

    assert all(layout in ['SINGLE', 'PAIRED'] for sample, layout in config['layout'].items())

    # if new samples were added, update the cache
    if len([sample for sample in samples.index if sample not in layout_cache]) is not 0:
        pickle.dump({**config['layout']}, open(layout_cachefile, "wb"))

# now only keep the layout of samples that are in samples.tsv
config['layout'] = {key: value for key, value in config['layout'].items() if key in samples.index}

logger.info("Done!\n\n")

# if samples are merged add the layout of the technical replicate to the config
if 'replicate' in samples and config.get('technical_replicates') == 'merge':
    for sample in samples.index:
        replicate = samples.loc[sample, 'replicate']
        config['layout'][replicate] = config['layout'][sample]

# after all is done, log (print) the configuration
logger.info("CONFIGURATION VARIABLES:")
for key, value in config.items():
     logger.info(f"{key: <23}: {value}")
logger.info("\n\n")

# save a copy of the latest samples and config file in the log_dir
# skip this step on Jenkins, as it runs in parallel
if os.getcwd() != config['log_dir'] and not os.getcwd().startswith('/var/lib/jenkins'):
    os.makedirs(config['log_dir'], exist_ok=True)
    for file in [config['samples']] + workflow.configfiles:
        src = os.path.join(os.getcwd(), file)
        dst = os.path.join(config['log_dir'], os.path.basename(file))
        shutil.copy(src, dst)


# check if a newer version of the Snakemake-workflows (master branch) is available
# if so, provide update instructions depending on the current branch
git_status = subprocess.check_output("""git fetch --dry-run -v 2>&1 | grep origin/master | cut -d "[" -f2 | cut -d "]" -f1""", shell=True).decode('ascii').strip()
if git_status != "up to date":
    current_branch = subprocess.check_output("""git branch | grep \* | awk '{print $2}'""", shell=True).decode('ascii').strip()
    cmd = " " if current_branch == "master" else " checkout master; git "
    logger.info(
        "A newer version of Snakemake-workflows is available!\n\n" +
        "To update, run:\n" +
        "\tgit" + cmd + "pull origin master\n\n"
    )
                      

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

# set global wildcard constraints (see workflow._wildcard_constraints)
wildcard_constraints:
    sample=any_given('sample'),
    sorting='coordinate|queryname',
    sorter='samtools|sambamba'

if 'assembly' in samples:
    wildcard_constraints:
        assembly=any_given('assembly'),

if 'replicate' in samples:
    wildcard_constraints:
        sample=any_given('sample', 'replicate'),
        replicate=any_given('replicate')             

if 'condition' in samples:
    wildcard_constraints:
        sample=any_given('sample', 'condition'),
        condition=any_given('condition')

if 'replicate' in samples and 'condition' in samples:
    wildcard_constraints:
        sample=any_given('sample', 'replicate', 'condition')


def get_workflow():
    return workflow.snakefile.split('/')[-2]


# set default parameters (parallel downloads and memory)
def convert_size(size_bytes, order=None):
    # https://stackoverflow.com/questions/5194057/better-way-to-convert-file-sizes-in-python/14822210#14822210
    size_name = ["b", "kb", "mb", "gb"]
    if order is None:
        order = int(math.floor(math.log(size_bytes, 1024)))
    s = int(size_bytes // math.pow(1024, order))
    return s, size_name[order]


# by default only one download in parallel (workflow fails on multiple on a single node)
workflow.global_resources = {**{'parallel_downloads': 3, 'deeptools_limit': 1, 'R_scripts': 1},
                             **workflow.global_resources}

# when the user specifies memory, use this and give a warning if it surpasses local memory
# (surpassing does not always have to be an issue -> cluster execution)
# if none specified set the memory to the max available on the computer
mem = psutil.virtual_memory()
if workflow.global_resources.get('mem_mb'):
    if workflow.global_resources['mem_mb'] > convert_size(mem.total, 2)[0]:
        logger.info(f"WARNING: The specified ram ({workflow.global_resources['mem_mb']} mb) surpasses the local machine\'s RAM ({convert_size(mem.total, 2)[0]} mb)")

    workflow.global_resources = {**workflow.global_resources,
                                 **{'mem_mb': np.clip(workflow.global_resources['mem_mb'], 0, convert_size(mem.total, 2)[0])}}
else:
    if workflow.global_resources.get('mem_gb', 0) > convert_size(mem.total, 3)[0]:
        logger.info(f"WARNING: The specified ram ({workflow.global_resources['mem_gb']} gb) surpasses the local machine\'s RAM ({convert_size(mem.total, 3)[0]} gb)")

    workflow.global_resources = {**workflow.global_resources,
                                 **{'mem_gb': np.clip(workflow.global_resources.get('mem_gb', 9999), 0, convert_size(mem.total, 3)[0])}}

def use_alignmentsieve(configdict):
    """
    helper function to check whether or not we use alignmentsieve
    """
    return configdict.get('min_mapping_quality', 0) > 0 or \
           configdict.get('tn5_shift', False) or \
           configdict.get('remove_blacklist', False) or \
           configdict.get('remove_mito', False)
