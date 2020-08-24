import contextlib
import genomepy
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
from functools import lru_cache

import norns
import numpy as np
import pandas as pd
import urllib.request
from bs4 import BeautifulSoup
from multiprocessing.pool import ThreadPool
from filelock import FileLock
from pandas_schema import Column, Schema
from pandas_schema.validation import MatchesPatternValidation, IsDistinctValidation

from snakemake.logging import logger
from snakemake.utils import validate
from snakemake.utils import min_version
from snakemake.exceptions import TerminatedException

import seq2science


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

assert all([col[0:7] not in ["Unnamed", ''] for col in samples]), \
    (f"\nEncountered unnamed column in {config['samples']}.\n" +
     f"Column names: {str(', '.join(samples.columns))}.\n")

# use pandasschema for checking if samples file is filed out correctly
allowed_pattern = r'[A-Za-z0-9_.\-%]+'
distinct_columns = ["sample"]
if "descriptive_name" in samples.columns:
    distinct_columns.append("descriptive_name")

distinct_schema = Schema(
    [Column(col, [MatchesPatternValidation(allowed_pattern),
                  IsDistinctValidation()] if col in distinct_columns else [MatchesPatternValidation(allowed_pattern)], allow_empty=True) for col in
     samples.columns])

errors = distinct_schema.validate(samples)

if len(errors):
    logger.error("\nThere are some issues with parsing the samples file:")
    for error in errors:
        logger.error(error)
    logger.error("")  # empty line
    raise TerminatedException

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
    if config.get('ignore_strandedness', True) or not any([field in list(samples['strandedness']) for field in ['yes', 'forward', 'reverse', 'no']]):
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


# check availability of assembly genomes and annotations


def get_workflow():
    return workflow.snakefile.split('/')[-2]


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


if "assembly" in samples:
    # control whether to custom extended assemblies
    if isinstance(config.get("custom_genome_extension"), str):
        config["custom_genome_extension"] = [config["custom_genome_extension"]]
    if isinstance(config.get("custom_annotation_extension"), str):
        config["custom_annotation_extension"] = [config["custom_annotation_extension"]]
    modified = config.get("custom_genome_extension") or config.get("custom_annotation_extension")
    all_assemblies = [assembly + "_custom" if modified else assembly for assembly in set(samples['assembly'])]
    suffix = "_custom" if modified else ""

    def list_providers(assembly):
        """
        Return a minimal list of providers to check
        """
        readme_file = os.path.join(config['genome_dir'], assembly, "README.txt")
        readme_provider = None
        if os.path.exists(readme_file):
            metadata, _ = genomepy.utils.read_readme(readme_file)
            readme_provider = metadata.get("provider", "").lower()

        if readme_provider in ["ensembl", "ucsc", "ncbi"]:
            providers = [readme_provider]
        elif config.get("provider"):
            providers = [config["provider"].lower()]
        else:
            providers = ["ensembl", "ucsc", "ncbi"]

        return providers


    def provider_with_file(file, assembly):
        """
        Returns the first provider which has the file for the assembly.
        file: annotation or genome
        """
        with open(os.devnull, "w") as null:
            with contextlib.redirect_stdout(null), contextlib.redirect_stderr(null):

                for provider in list_providers(assembly):
                    p = genomepy.ProviderBase.create(provider)
                    if assembly in p.genomes:
                        if (file == "annotation" and p.get_annotation_download_link(assembly)) \
                                or (file == "genome" and p.get_genome_download_link(assembly)):
                            return provider
        return None


    # determine provider for each new assembly
    providersfile = os.path.expanduser('~/.config/seq2science/providers.p')
    providersfile_lock = os.path.expanduser('~/.config/seq2science/providers.p.lock')
    prep_filelock(providersfile_lock, 30)
    with FileLock(providersfile_lock):
        providers = dict()
        if os.path.exists(providersfile):
            providers = pickle.load(open(providersfile, "rb"))

        if any([assembly not in providers for assembly in set(samples["assembly"])]):
            logger.info("Determining assembly providers")

            for assembly in set(samples["assembly"]):
                if assembly not in providers:
                    file = os.path.join(config['genome_dir'], assembly, assembly)
                    providers[assembly] = {"genome": None, "annotation": None}

                    # check if genome and annotations exist locally
                    if os.path.exists(f"{file}.fa"):
                        providers[assembly]["genome"] = "local"
                    if all(os.path.exists(f) for f in [f"{file}.annotation.gtf", f"{file}.annotation.bed"]):
                        providers[assembly]["annotation"] = "local"

                    # check if the annotation can be downloaded
                    if providers[assembly]["annotation"] is None:
                        annotion_provider = provider_with_file("annotation", assembly)
                        if annotion_provider:
                            providers[assembly]["genome"] = annotion_provider  # exists if annotation does
                            providers[assembly]["annotation"] = annotion_provider

                    # check if the genome can be downloaded
                    if providers[assembly]["genome"] is None:
                        genome_provider = provider_with_file("genome", assembly)
                        providers[assembly]["genome"] = genome_provider

            pickle.dump(providers, open(providersfile, "wb"))

    # check the providers for the required assemblies
    annotation_required = "rna_seq" in get_workflow() or config["aligner"] == "star"
    for assembly in set(samples["assembly"]):
        file = os.path.join(config['genome_dir'], assembly, assembly)
        if providers[assembly]["genome"] is None and not os.path.exists(f"{file}.fa"):
            logger.info(
                f"Could not download assembly {assembly}.\n"
                f"Find alternative assemblies with `genomepy search {assembly}`"
            )
            exit(1)

        if providers[assembly]["annotation"] is None and \
                not all(os.path.exists(f) for f in [f"{file}.annotation.gtf", f"{file}.annotation.bed"]):
            logger.info(
                f"No annotation for assembly {assembly} can be downloaded. Another provider (and "
                f"thus another assembly name) might have gene annotations.\n"
                f"Find alternative assemblies with `genomepy search {assembly}`"
            )
            if annotation_required:
                exit(1)
            time.sleep(0 if config.get("debug") else 2)  # give some time to read the message


    def ori_assembly(assembly):
        """
        remove _SI from assembly if is was added
        """
        return assembly[:-7] if assembly.endswith("_custom") and modified else assembly


    @lru_cache(maxsize=None)
    def has_annotation(assembly):
        """
        Returns True/False on whether or not the assembly has an annotation.
        """
        return True if providers[ori_assembly(assembly)]["annotation"] else False


# sample layouts


# check if a sample is single-end or paired end, and store it
logger.info("Checking if samples are single-end or paired-end...")
logger.info("This can take some time.")
layout_cachefile = os.path.expanduser('~/.config/seq2science/layouts.p')
layout_cachefile_lock = os.path.expanduser('~/.config/seq2science/layouts.p.lock')
sample_to_ena_single_url = os.path.expanduser('~/.config/seq2science/ena_single.p')
sample_to_ena_paired_url = os.path.expanduser('~/.config/seq2science/ena_paired.p')
sample_to_ena_url_lock = os.path.expanduser('~/.config/seq2science/ena.p.lock')


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


def get_layout_trace1(sample, timeout=10, max_tries=1):
    """
    Parse the ncbi trace website to check if a read has 1, 2, or 3 spots.
    Will fail if sample is not on ncbi database, however does not have the problem that uploader
    filled out the form wrongly.
    Complementary method of get_layout_eutils
    """
    for i in range(max_tries):
        try:
            url = f"https://www.ncbi.nlm.nih.gov/sra/?term={sample}"

            conn = urllib.request.urlopen(url, timeout=timeout)
            html = conn.read()

            soup = BeautifulSoup(html, features="html.parser")
            links = soup.find_all('a')

            vals = []
            for tag in links:
                link = tag.get('href', None)
                if link is not None and 'SRR' in link:
                    SRR = link[link.find("SRR"):]
                    trace_conn = urllib.request.urlopen("https:" + link, timeout=timeout)
                    trace_html = trace_conn.read()
                    x = re.search("This run has (\d) read", str(trace_html))

                    # if there are spots without info, then just ignore this sample
                    if len(re.findall(", average length: 0", str(trace_html))) > 0:
                        break
                    elif x.group(1) == '1':
                        vals.append(('SINGLE', SRR))
                    elif x.group(1) == '2':
                        vals.append(('PAIRED', SRR))
            if len(vals) > 0:
                return vals
        except:
            pass
    return None


def get_layout_trace2(sample, timeout=10, max_tries=1):
    """
    Yet another sample lookup fallback.
    """
    for i in range(max_tries):
        try:
            url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sample}"

            conn = urllib.request.urlopen(url, timeout=timeout)
            html = conn.read()

            soup = BeautifulSoup(html, features="html.parser")
            links = soup.find_all('a')

            for tag in links:
                link = tag.get('href', None)
                if link is not None and 'SRX' in link:
                    SRX = link[link.find("SRX"):]
                    return get_layout_trace1(SRX, timeout)
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

    trace_layout1 = {}
    trace_layout2 = {}
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
            config['layout'][sample] = 'SINGLE'
        elif all(os.path.exists(path) for path in expand(f'{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
            config['layout'][sample] = 'PAIRED'
        elif sample.startswith(('GSM', 'SRR', 'ERR', 'DRR')):
            eutils_layout[sample] = eutils_tp.apply_async(get_layout_eutils, (sample,))
            trace_layout1[sample] = trace_tp.apply_async(get_layout_trace1, (sample, 10, 3))
            trace_layout2[sample] = trace_tp.apply_async(get_layout_trace2, (sample, 10, 3))

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
                        **{k: v.get()[0][0] for k, v in trace_layout1.items() if v.get() is not None},
                        **{k: v.get()[0][0] for k, v in trace_layout2.items() if v.get() is not None}}

    assert all(layout in ['SINGLE', 'PAIRED'] for sample, layout in config['layout'].items())

    # if new samples were added, update the cache
    if len([sample for sample in samples.index if sample not in layout_cache]) is not 0:
        pickle.dump({**config['layout']}, open(layout_cachefile, "wb"))


# now only keep the layout of samples that are in samples.tsv
config['layout'] = {**{key: value for key, value in config['layout'].items() if key in samples.index},
                    **{key: value for key, value in config['layout'].items() if "control" in samples and key in samples["control"].values}}

bad_samples = [sample for sample in samples.index if sample not in config["layout"]]
if len(bad_samples) > 0:
    logger.error(f"\nThe instructions to lookup sample(s) {' '.join(bad_samples)} online failed!\n"
                 f"Are you sure these sample(s) exists..? Downloading samples with restricted "
                 f"access is currently not supported. We advise you to download the sample "
                 f"manually, and continue the pipeline from there on.\n")
    raise TerminatedException

sample_to_srr = {**{k: v.get() for k, v in trace_layout1.items() if v.get() is not None},
                 **{k: v.get() for k, v in trace_layout2.items() if v.get() is not None}}

trace_tp.close()
eutils_tp.close()

def url_is_alive(url):
    """
    Checks that a given URL is reachable.
    https://gist.github.com/dehowell/884204
    """
    for i in range(3):
        try:
            request = urllib.request.Request(url)
            request.get_method = lambda: 'HEAD'

            urllib.request.urlopen(request, timeout=5)
            return True
        except:
            continue
    return False

logger.info("Done!\n\n")
logger.info("Now checking if the samples are on the ENA database..")
logger.info("This can also take some time!")

ena_single_end_urls = dict()
ena_paired_end_urls = dict()

# trace_tp = ThreadPool(40)

# now check if we can simply download the fastq from ENA
prep_filelock(sample_to_ena_url_lock, 5*60)

with FileLock(sample_to_ena_url_lock):
    # try to load the layout cache, otherwise defaults to empty dictionary
    try:
        ena_single_end_urls = pickle.load(open(sample_to_ena_single_url, "rb"))
        ena_paired_end_urls = pickle.load(open(sample_to_ena_paired_url, "rb"))
    except FileNotFoundError:
        ena_single_end_urls = {}
        ena_paired_end_urls = {}

    for sample in samples.index:
        # do not check if in cache
        if sample in ena_single_end_urls or sample in ena_paired_end_urls:
            continue

        # do not check if the file already exists
        if os.path.exists(expand(f'{{fastq_dir}}/{sample}.{{fqsuffix}}.gz', **config)[0]) or \
           all(os.path.exists(path) for path in expand(f'{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
            continue

        srrs = sample_to_srr.get(sample, None)
        if srrs is not None:
            layout = [srr[0] for srr in srrs]
            if len(set(layout)) > 1:
                raise ValueError("I can not deal with mixes of samples!")
            if len(set(layout)) == 0:
                continue
            layout = layout[0]

            for srr in [srr[1] for srr in srrs]:
                prefix = srr[:6]
                suffix = f"/{int(srr[9:]):03}" if len(srr) >= 10 else ""

                fasp_address = "era-fasp@fasp.sra.ebi.ac.uk:"
                wget_address = "ftp://ftp.sra.ebi.ac.uk/"

                if layout == "SINGLE":
                    wget_url = f"{wget_address}vol1/fastq/{prefix}{suffix}/{srr}/{srr}.fastq.gz"
                    fasp_url = f"{fasp_address}vol1/fastq/{prefix}{suffix}/{srr}/{srr}.fastq.gz"
                    if url_is_alive(wget_url):
                        url = fasp_url if config.get("ascp_path") else wget_url
                        ena_single_end_urls.setdefault(sample, []).append((srr, url))
                elif layout == "PAIRED":
                    wget_urls = [f"{wget_address}vol1/fastq/{prefix}{suffix}/{srr}/{srr}_1.fastq.gz",
                                 f"{wget_address}vol1/fastq/{prefix}{suffix}/{srr}/{srr}_2.fastq.gz"]
                    fasp_urls = [f"{fasp_address}vol1/fastq/{prefix}{suffix}/{srr}/{srr}_1.fastq.gz",
                                 f"{fasp_address}vol1/fastq/{prefix}{suffix}/{srr}/{srr}_2.fastq.gz"]
                    if all(url_is_alive(url) for url in wget_urls):
                        urls = fasp_urls if config.get("ascp_path") else wget_urls
                        ena_paired_end_urls.setdefault(sample, []).append((srr, urls))
                else:
                    raise NotImplementedError

    pickle.dump(ena_single_end_urls, open(sample_to_ena_single_url, "wb"))
    pickle.dump(ena_paired_end_urls, open(sample_to_ena_paired_url, "wb"))


logger.info("Done!\n\n")

# if samples are merged add the layout of the technical replicate to the config
if 'replicate' in samples:
    for sample in samples.index:
        replicate = samples.loc[sample, 'replicate']
        config['layout'][replicate] = config['layout'][sample]


# workflow


def any_given(*args, prefix="", suffix=""):
    """
    returns a regex compatible string of all elements in the samples.tsv column given by the input
    """
    elements = []
    for column_name in args:
        if column_name is 'sample':
            elements.extend(samples.index)
        elif column_name in samples:
            elements.extend(samples[column_name])

    elements = [prefix + element + suffix for element in elements if isinstance(element, str)]
    return '|'.join(set(elements))

# set global wildcard constraints (see workflow._wildcard_constraints)
sample_constraints = ["sample"]
wildcard_constraints:
    sorting='coordinate|queryname',
    sorter='samtools|sambamba',

if 'assembly' in samples:
    wildcard_constraints:
        raw_assembly=any_given('assembly'),
        assembly=any_given('assembly', suffix=config["spike_suffix"] if modified else ""),

if 'replicate' in samples:
    sample_constraints.append("replicate")
    wildcard_constraints:
        replicate=any_given('replicate')

if 'condition' in samples:
    sample_constraints.append("condition")
    wildcard_constraints:
        condition=any_given('condition')

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
