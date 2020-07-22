import math
import os.path
import psutil
import pickle
import re
import shutil
import subprocess
import time
import copy
import json
import requests

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


logger.info(
"""\
               ____  ____   __              
              / ___)(  __) /  \             
              \___ \ ) _) (  O )            
              (____/(____) \__\)            
                     ____                   
                    (___ \                  
                     / __/                  
                    (____)                  
   ____   ___   __   ____  __ _   ___  ____ 
  / ___) / __) (  ) (  __)(  ( \ / __)(  __)
  \___ \( (__   )(   ) _) /    /( (__  ) _) 
  (____/ \___) (__) (____)\_)__) \___)(____)

docs: https://vanheeringen-lab.github.io/seq2science
"""
)

include: f"{config['rule_dir']}/explain.smk"


if workflow.conda_frontend == "conda":
    logger.info("NOTE: seq2science is using the conda frontend, for faster environment creation install mamba.")
# give people a second to appreciate this beautiful ascii art
time.sleep(1)


# config.yaml(s)


# check config file for correct directory names
for key, value in config.items():
    if '_dir' in key:
        assert value not in [None, '', ' '], f"\n{key} cannot be empty. For current directory, set '{key}: .'\n"
        # allow tilde as the first character, \w = all letters and numbers
        # ignore on Jenkins. it puts @ in the path
        assert (re.match('^[~\w_./-]*$', value[0]) and re.match('^[\w_./-]*$', value[1:])) \
               or os.getcwd().startswith('/var/lib/jenkins'), \
            (f"\nIn the config.yaml you set '{key}' to '{value}'. " +
              "Please use file paths that only contain letters, " +
              "numbers and any of the following symbols: underscores (_), periods (.), " +
              "and minuses (-). The first character may also be a tilde (~).\n")
        config[key] = os.path.expanduser(value)

# make sure that difficult data-types (yaml objects) are in correct data format
for kw in ['aligner', 'quantifier', 'bam_sorter']:
    if isinstance(config.get(kw, None), str):
        config[kw] = {config[kw]: {}}

# validate and complement the config dict
for schema in config_schemas:
    validate(config, schema=f"{config['rule_dir']}/../schemas/config/{schema}.schema.yaml")

# check if paired-end filename suffixes are lexicographically ordered
config['fqext'] = [config['fqext1'], config['fqext2']]
assert sorted(config['fqext'])[0] == config['fqext1'], \
    ("\nThe paired-end filename suffixes must be lexicographically ordered!\n" +
     f"Example suffixes: fqext1: R1, fqext2: R2\n" +
     f"Your suffixes:    fqext1: {config['fqext1']}, fqext2: {config['fqext2']}\n")

# read the config.yaml (not the profile)
user_config = norns.config(config_file=workflow.overwrite_configfiles[0])

# make absolute paths, nest default dirs in result_dir and cut off trailing slashes
config['result_dir'] = re.split("\/$", os.path.abspath(config['result_dir']))[0]
config['samples'] = os.path.abspath(config['samples'])
for key, value in config.items():
    if key.endswith("_dir"):
        if key in ['result_dir', 'genome_dir', 'rule_dir'] or key in user_config:
            value = os.path.abspath(value)
        else:
            value = os.path.abspath(os.path.join(config['result_dir'], value))
        config[key] = re.split("\/$", value)[0]


# samples.tsv


# read the samples.tsv file as all text, drop comment lines
samples = pd.read_csv(config["samples"], sep='\t', dtype='str', comment='#')
samples.columns = samples.columns.str.strip()

# check that the columns are named
assert all([col[0:7] not in ["Unnamed", ''] for col in samples]), \
    (f"\nEncountered unnamed column in {config['samples']}.\n" +
     f"Column names: {str(', '.join(samples.columns))}.\n")

# check that the columns contains no irregular characters
assert not any(samples.columns.str.contains('[^A-Za-z0-9_.\-%]+', regex=True)), \
    (f"\n{config['samples']} may only contain letters, numbers and " +
    "percentage signs (%), underscores (_), periods (.), or minuses (-).\n")

# check that the file contains no irregular characters
assert not any([any(samples[col].str.contains('[^A-Za-z0-9_.\-%]+', regex=True, na=False)) for col in samples if col != "control"]), \
    (f"\n{config['samples']} may only contain letters, numbers and " +
    "percentage signs (%), underscores (_), periods (.), or minuses (-).\n")

# check that sample names are unique
assert len(samples["sample"]) == len(set(samples["sample"])), \
    (f"\nDuplicate samples found in {config['samples']}:\n" +
     f"{samples[samples.duplicated(['sample'], keep=False)].to_string()}\n")

# for each column, if found in samples.tsv:
# 1) if it is incomplete, fill the blanks with replicate/sample names
# (sample names if replicates are not found/applicable)
# 2) drop column if it is identical to the replicate/sample column, or if not needed
if 'replicate' in samples:
    samples['replicate'] = samples['replicate'].mask(pd.isnull, samples['sample'])
    if samples['replicate'].tolist() == samples['sample'].tolist() or config.get('technical_replicates') == 'keep':
        samples = samples.drop(columns=['replicate'])
if 'condition' in samples:
    samples['condition'] = samples['condition'].mask(pd.isnull, samples['replicate']) if 'replicate' in samples else \
        samples['condition'].mask(pd.isnull, samples['sample'])
    if samples['condition'].tolist() == samples['sample'].tolist() or config.get('biological_replicates') == 'keep':
        samples = samples.drop(columns=['condition'])
if 'descriptive_name' in samples:
    samples['descriptive_name'] = samples['descriptive_name'].mask(pd.isnull, samples['replicate']) if \
        'replicate' in samples else samples['descriptive_name'].mask(pd.isnull, samples['sample'])
    if ('replicate' in samples and samples['descriptive_name'].to_list() == samples['replicate'].to_list()) or \
        samples['descriptive_name'].to_list() == samples['sample'].to_list():
        samples = samples.drop(columns=['descriptive_name'])
if 'strandedness' in samples:
    samples['strandedness'] = samples['strandedness'].mask(pd.isnull, 'nan')
    if config['filter_bam_by_strand'] is False or not any([field in list(samples['strandedness']) for field in ['forward', 'yes', 'reverse']]):
        samples = samples.drop(columns=['strandedness'])

if 'replicate' in samples:
    # check if replicate names are unique between assemblies
    r = samples[['assembly', 'replicate']].drop_duplicates().set_index('replicate')
    for replicate in r.index:
        assert len(r[r.index == replicate]) == 1, \
            ("\nReplicate names must be different between assemblies.\n" +
             f"Replicate name '{replicate}' was found in assemblies {r[r.index == replicate]['assembly'].tolist()}.")

# check if sample, replicate and condition names are unique between the columns
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
    validate(samples, schema=f"{config['rule_dir']}/../schemas/samples/{schema}.schema.yaml")

sanitized_samples = copy.copy(samples)

samples = samples.set_index('sample')
samples.index = samples.index.map(str)


# sample layouts
# check if a sample is single-end or paired end, and store it
logger.info("Checking if samples are single-end or paired-end...")
layout_cachefile = os.path.expanduser('~/.config/seq2science/layouts.p')
layout_cachefile_lock = os.path.expanduser('~/.config/seq2science/layouts.p.lock')

def prep_filelock(lock_file, max_age=10):
    """
    create the directory for the lock_file if needed
    and remove locks older than the max_age (in seconds)
    """
    os.makedirs(os.path.dirname(lock_file), exist_ok=True)

    # sometimes two jobs start in parallel and try to delete at the same time
    try:
        # ignore locks that are older than the max_age
        if os.path.exists(lock_file) and \
                time.time() - os.stat(lock_file).st_mtime > max_age:
            os.unlink(lock_file)
    except FileNotFoundError:
         pass


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
        return None


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
prep_filelock(layout_cachefile_lock, 5*60)

with FileLock(layout_cachefile_lock):
    # try to load the layout cache, otherwise defaults to empty dictionary
    try:
        layout_cache = pickle.load(open(layout_cachefile, "rb"))
    except FileNotFoundError:
        layout_cache = {}


    trace_tp = ThreadPool(20)
    eutils_tp = ThreadPool(config.get('ncbi_requests', 3) // 2)

    trace_layout = {}
    eutils_layout = {}
    config['layout'] = {}

    # now do a request for each sample that was not in the cache
    all_samples = [sample for sample in samples.index if sample not in layout_cache]
    if "control" in samples:
        for control in set(samples["control"]):
            if control not in layout_cache and isinstance(control, str):  # ignore nans
                all_samples.append(control)

    for sample in all_samples:
        if os.path.exists(expand(f'{{fastq_dir}}/{sample}.{{fqsuffix}}.gz', **config)[0]):
            config['layout'][sample] ='SINGLE'
        elif all(os.path.exists(path) for path in expand(f'{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
            config['layout'][sample] ='PAIRED'
        elif sample.startswith(('GSM', 'SRR', 'ERR', 'DRR')):
            eutils_layout[sample] = eutils_tp.apply_async(get_layout_eutils, (sample,))
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
                        **{k: v for k, v in config['layout'].items()},
                        **{k: v.get() for k, v in eutils_layout.items() if v.get() is not None},
                        **{k: v.get() for k, v in trace_layout.items() if v.get() is not None}}

    assert all(layout in ['SINGLE', 'PAIRED'] for sample, layout in config['layout'].items())

    # if new samples were added, update the cache
    if len([sample for sample in samples.index if sample not in layout_cache]) is not 0:
        pickle.dump({**config['layout']}, open(layout_cachefile, "wb"))


# now only keep the layout of samples that are in samples.tsv
config['layout'] = {**{key: value for key, value in config['layout'].items() if key in samples.index},
                    **{key: value for key, value in config['layout'].items() if "control" in samples and key in samples["control"].values}}

for sample in samples.index:
    if sample not in config["layout"]:
        raise ValueError(f"The command to lookup sample {sample} online failed!\n"
                         f"Are you sure this sample exists..? Downloading samples with restricted "
                         f"access is currently not supported. We advise you to download the sample "
                         f"manually, and continue the pipeline from there on.")


logger.info("Done!\n\n")

# if samples are merged add the layout of the technical replicate to the config
if 'replicate' in samples:
    for sample in samples.index:
        replicate = samples.loc[sample, 'replicate']
        config['layout'][replicate] = config['layout'][sample]


# workflow


def get_workflow():
    return workflow.snakefile.split('/')[-2]


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

    elements = [element for element in elements if isinstance(element, str)]
    return '|'.join(set(elements))

# set global wildcard constraints (see workflow._wildcard_constraints)
sample_constraints = ["sample"]
wildcard_constraints:
    sorting='coordinate|queryname',
    sorter='samtools|sambamba'

if 'assembly' in samples:
    wildcard_constraints:
        assembly=any_given('assembly'),

if 'replicate' in samples:
    sample_constraints = ["sample", "replicate"]
    wildcard_constraints:
        replicate=any_given('replicate')

if 'condition' in samples:
    sample_constraints = ["sample", "condition"]
    wildcard_constraints:
        condition=any_given('condition')

if 'replicate' in samples and 'condition' in samples:
    sample_constraints = ["sample", "replicate", "condition"]

if "control" in samples:
    sample_constraints.append("control")

wildcard_constraints:
    sample=any_given(*sample_constraints)


# set default parameters (parallel downloads and memory)
def convert_size(size_bytes, order=None):
    # https://stackoverflow.com/questions/5194057/better-way-to-convert-file-sizes-in-python/14822210#14822210
    size_name = ["b", "kb", "mb", "gb"]
    if order is None:
        order = int(math.floor(math.log(size_bytes, 1024)))
    s = int(size_bytes // math.pow(1024, order))
    return s, size_name[order]

# make sure the snakemake version corresponds to version in environment
min_version("5.18")

# set some defaults
workflow.global_resources = {**{'parallel_downloads': 3, 'deeptools_limit': 16, 'R_scripts': 1},
                             **workflow.global_resources}

# when the user specifies memory, use this and give a warning if it surpasses local memory
# (surpassing does not always have to be an issue -> cluster execution)
# if none specified set the memory to the max available on the computer
mem = psutil.virtual_memory()
if workflow.global_resources.get('mem_mb'):
    if workflow.global_resources['mem_mb'] > convert_size(mem.total, 2)[0]:
        logger.info(f"WARNING: The specified ram ({workflow.global_resources['mem_mb']} mb) surpasses the local machine\'s RAM ({convert_size(mem.total, 2)[0]} mb)")

    workflow.global_resources = {**{'mem_mb': np.clip(workflow.global_resources['mem_mb'], 0, convert_size(mem.total, 2)[0])},
                                 **workflow.global_resources}
else:
    if workflow.global_resources.get('mem_gb', 0) > convert_size(mem.total, 3)[0]:
        logger.info(f"WARNING: The specified ram ({workflow.global_resources['mem_gb']} gb) surpasses the local machine\'s RAM ({convert_size(mem.total, 3)[0]} gb)")

    workflow.global_resources = {**{'mem_gb': np.clip(workflow.global_resources.get('mem_gb', 9999), 0, convert_size(mem.total, 3)[0])},
                                 **workflow.global_resources}

# record which assembly trackhubs are found on UCSC
if config.get("create_trackhub"):
    hubfile = os.path.expanduser('~/.config/seq2science/ucsc_trackhubs.p')
    hubfile_lock = os.path.expanduser('~/.config/seq2science/ucsc_trackhubs.p.lock')
    prep_filelock(hubfile_lock)

    with FileLock(hubfile_lock):
        if not os.path.exists(hubfile):
            # check for response of ucsc
            response = requests.get(f"https://genome.ucsc.edu/cgi-bin/hgGateway",
                                    allow_redirects=True)
            assert response.ok, "Make sure you are connected to the internet"

            with urllib.request.urlopen("https://api.genome.ucsc.edu/list/ucscGenomes") as url:
                data = json.loads(url.read().decode())['ucscGenomes']

            # generate a dict ucsc assemblies
            ucsc_assemblies = dict()
            for key, values in data.items():
                ucsc_assemblies[key.lower()] = [key, values.get("description", "")]

            # save to file
            pickle.dump(ucsc_assemblies, open(hubfile, "wb"))

        # read hubfile
        ucsc_assemblies = pickle.load(open(hubfile, "rb"))

onstart:
    # save a copy of the latest samples and config file(s) in the log_dir
    # skip this step on Jenkins, as it runs in parallel
    if os.getcwd() != config['log_dir'] and not os.getcwd().startswith('/var/lib/jenkins'):
        os.makedirs(config['log_dir'], exist_ok=True)
        for n, file in enumerate([config['samples']] + workflow.overwrite_configfiles):
            src = os.path.join(os.getcwd(), file)
            dst = os.path.join(config['log_dir'], os.path.basename(file) if n<2 else "profile.yaml")
            shutil.copy(src, dst)
onsuccess:
    if config.get("email") not in ["none@provided.com", "yourmail@here.com", None]:
        os.system(f"""echo "Succesful pipeline run! :)" | mail -s "The seq2science pipeline finished succesfully." {config["email"]} 2> /dev/null""")
onerror:
    if config.get("email") not in ["none@provided.com", "yourmail@here.com", None]:
        os.system(f"""echo "Unsuccessful pipeline run! :(" | mail -s "The seq2science pipeline finished prematurely..." {config["email"]} 2> /dev/null """)

include: "../rules/configuration_workflows.smk"
