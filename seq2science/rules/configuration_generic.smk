import contextlib
import genomepy
import math
import os.path
import pickle
import re
import time
import copy
import json
import requests
from functools import lru_cache
import yaml

import xdg
import numpy as np
import pandas as pd
import urllib.request
from filelock import FileLock
from pandas_schema import Column, Schema
from pandas_schema.validation import MatchesPatternValidation, IsDistinctValidation

from snakemake.logging import logger
from snakemake.utils import validate
from snakemake.utils import min_version
from snakemake.exceptions import TerminatedException

import seq2science
from seq2science.util import samples2metadata, prep_filelock, url_is_alive, color_parser



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


if workflow.conda_frontend == "conda":
    logger.info("NOTE: seq2science is using the conda frontend, for faster environment creation install mamba.")
# give people a second to appreciate this beautiful ascii art
time.sleep(1)

# get the cache and config dirs
CACHE_DIR = os.path.join(xdg.XDG_CACHE_HOME, "seq2science", seq2science.__version__)

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
for kw in ['aligner', 'quantifier', 'bam_sorter', "trimmer"]:
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
with open(workflow.overwrite_configfiles[0], 'r') as stream:
    user_config = yaml.safe_load(stream)

# make absolute paths, nest default dirs in result_dir and cut off trailing slashes
config['result_dir'] = re.split("\/$", os.path.abspath(config['result_dir']))[0]

if not url_is_alive(config['samples']):
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
if 'colors' in samples:
    samples['colors'] = samples['colors'].mask(pd.isnull, '0,0,0')  # nan -> black
    samples['colors'] = [color_parser(c) for c in samples['colors']]  # convert input to HSV color
    if not config.get('create_trackhub', False):
        samples = samples.drop(columns=['colors'])

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

sequencing_protocol = get_workflow()\
    .replace('alignment',  'Alignment')\
    .replace('atac_seq',   'ATAC-seq')\
    .replace('chip_seq',   'ChIP-seq')\
    .replace('rna_seq',    'RNA-seq')\
    .replace('scatac_seq', 'scATAC-seq')


if "assembly" in samples:
    # control whether to custom extended assemblies
    if isinstance(config.get("custom_genome_extension"), str):
        config["custom_genome_extension"] = [config["custom_genome_extension"]]
    if isinstance(config.get("custom_annotation_extension"), str):
        config["custom_annotation_extension"] = [config["custom_annotation_extension"]]
    modified = config.get("custom_genome_extension") or config.get("custom_annotation_extension")
    all_assemblies = [assembly + config["custom_assembly_suffix"] if modified else assembly for assembly in set(samples['assembly'])]
    suffix = config["custom_assembly_suffix"] if modified else ""

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
    providersfile = f"{CACHE_DIR}/providers.p"
    providersfile_lock = f"{CACHE_DIR}/providers.p.lock"
    for _ in range(2):
        # we get two tries, in case parallel executions are interfering with one another
        try:
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
                            if any(os.path.exists(f) for f in [f"{file}.annotation.gtf", f"{file}.annotation.gtf.gz"]) and \
                                any(os.path.exists(f) for f in [f"{file}.annotation.bed", f"{file}.annotation.bed.gz"]):
                                providers[assembly]["annotation"] = "local"

                            # check if the annotation can be downloaded
                            if providers[assembly]["annotation"] is None:
                                annotion_provider = provider_with_file("annotation", assembly)
                                if annotion_provider:
                                    providers[assembly]["genome"] = annotion_provider  # genome always exists if annotation does
                                    providers[assembly]["annotation"] = annotion_provider

                            # check if the genome can be downloaded
                            if providers[assembly]["genome"] is None:
                                genome_provider = provider_with_file("genome", assembly)
                                providers[assembly]["genome"] = genome_provider

                    pickle.dump(providers, open(providersfile, "wb"))
                break
        except FileNotFoundError:
            time.sleep(1)
    else:
        logger.error("There were some problems with locking the seq2science cache. Please try again in a bit.")
        raise TerminatedException

    # check the providers for the required assemblies
    annotation_required = "rna_seq" in get_workflow() or config["aligner"] == "star"
    _has_annot = dict()
    for assembly in set(samples["assembly"]):
        file = os.path.join(config['genome_dir'], assembly, assembly)
        if providers[assembly]["genome"] is None and not os.path.exists(f"{file}.fa"):
            logger.info(
                f"Could not download assembly {assembly}.\n"
                f"Find alternative assemblies with `genomepy search {assembly}`"
            )
            exit(1)

        if providers[assembly]["annotation"] is None and not (
                any(os.path.exists(f) for f in [f"{file}.annotation.gtf", f"{file}.annotation.gtf.gz"]) and
                any(os.path.exists(f) for f in [f"{file}.annotation.bed", f"{file}.annotation.bed.gz"])
        ):
            logger.info(
                f"No annotation for assembly {assembly} can be downloaded. Another provider (and "
                f"thus another assembly name) might have a gene annotation.\n"
                f"Find alternative assemblies with `genomepy search {assembly}`\n"
            )
            if annotation_required:
                exit(1)
            _has_annot[assembly] = False
            time.sleep(0 if config.get("debug") else 2)  # give some time to read the message
        else:
            _has_annot[assembly] = True


    def ori_assembly(assembly):
        """
        remove the extension suffix from an assembly if is was added.
        """
        return assembly[:-len(config["custom_assembly_suffix"])] if \
            assembly.endswith(config["custom_assembly_suffix"]) and modified else assembly


    def custom_assembly(assembly):
        """
        add extension suffix to an assembly if is wasn't yet added.
        """
        return assembly if assembly.endswith(config["custom_assembly_suffix"]) else \
         (assembly + config["custom_assembly_suffix"])


    @lru_cache(maxsize=None)
    def has_annotation(assembly):
        """
        Returns True/False on whether or not the assembly has an annotation.
        """
        return _has_annot[assembly]

else:
    modified = False


# sample layouts


# check if a sample is single-end or paired end, and store it
logger.info("Checking if samples are available online...")
logger.info("This can take some time.")

# make a collection of all samples
all_samples = [sample for sample in samples.index]
if "control" in samples:
    for control in set(samples["control"]):
        if isinstance(control, str):  # ignore nans
            all_samples.append(control)

pysradb_cache = f"{CACHE_DIR}/pysradb.p"
pysradb_cache_lock = f"{CACHE_DIR}/pysradb.p.lock"
for _ in range(2):
    # we get two tries, in case parallel executions are interfering with one another
    try:
        prep_filelock(pysradb_cache_lock, 30)
        with FileLock(pysradb_cache_lock):
            try:
                sampledict = pickle.load(open(pysradb_cache, "rb"))
            except FileNotFoundError:
                sampledict = {}

            missing_samples = [sample for sample in all_samples if sample not in sampledict.keys()]
            if len(missing_samples) > 0:
                sampledict.update(samples2metadata(missing_samples, config, logger))

            pickle.dump(sampledict, open(pysradb_cache, "wb"))

            # only keep samples for this run
            sampledict = {sample: values for sample, values in sampledict.items() if sample in all_samples}
        break
    except FileNotFoundError:
        time.sleep(1)
else:
    logger.error("There were some problems with locking the seq2science cache. Please try again in a bit.")
    raise TerminatedException

logger.info("Done!\n\n")

# now check where to download which sample
ena_single_end = [run for values in sampledict.values() if (values["layout"] == "SINGLE") and values.get("ena_fastq_ftp") is not None for run in values["runs"]]
ena_paired_end = [run for values in sampledict.values() if (values["layout"] == "PAIRED") and values.get("ena_fastq_ftp") is not None for run in values["runs"]]

# get download link per run
run2download = dict()
for sample, values in sampledict.items():
    for run in values.get("runs", []):
        if values["ena_fastq_ftp"] and values["ena_fastq_ftp"][run]:
            if not (config.get("ascp_path") and config.get("ascp_key")):
                run2download[run] = [url.replace("era-fasp@fasp", "ftp") for url in values["ena_fastq_ftp"][run]]
            else:
                run2download[run] = values["ena_fastq_ftp"][run]

# if samples are merged add the layout of the technical replicate to the config
failed_samples = dict()
if 'replicate' in samples:
    for sample in samples.index:
        replicate = samples.loc[sample, 'replicate']
        if replicate not in sampledict:
            sampledict[replicate] = {'layout':  sampledict[sample]['layout']}
        elif sampledict[replicate]['layout'] != sampledict[sample]['layout']:
            assembly = samples.loc[sample, "assembly"]
            treps = samples[(samples["assembly"] == assembly) & (samples["replicate"] == replicate)].index
            failed_samples.setdefault(replicate, set()).update({trep for trep in treps}) 

if len(failed_samples):
    logger.error("Your technical replicates consist of a mix of single-end and paired-end samples!")
    logger.error("This is not supported.\n")

    for replicate, samples in failed_samples.items():
        logger.error(f"{replicate}:")
        for sample in samples:
            logger.error(f"\t{sample}: {sampledict[sample]['layout']}")
        logger.error("\n")
    raise TerminatedException

# workflow


def any_given(*args, prefix="", suffix=""):
    """
    returns a regex compatible string of all elements in the samples.tsv column given by the input
    """
    elements = []
    for column_name in args:
        if column_name == 'sample':
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
        assembly=any_given('assembly', suffix=config["custom_assembly_suffix"] if modified else ""),

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


# make sure the snakemake version corresponds to version in environment
min_version("5.18")


# record which assembly trackhubs are found on UCSC
if config.get("create_trackhub"):
    hubfile = f"{CACHE_DIR}/ucsc_trackhubs.p"
    hubfile_lock = f"{CACHE_DIR}/ucsc_trackhubs.p.lock"
    for _ in range(2):
        # we get two tries, in case parallel executions are interfering with one another
        try:
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
            break
        except FileNotFoundError:
            time.sleep(1)
    else:
        logger.error("There were some problems with locking the seq2science cache. Please try again in a bit.")
        raise TerminatedException

