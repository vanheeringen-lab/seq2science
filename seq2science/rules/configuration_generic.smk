"""
all logic not related to any specific workflows should be here.
"""

import os.path
import pickle
import re
import time
import copy
import json
import requests
from functools import lru_cache
import yaml
import sys

import xdg
import pandas as pd
import urllib.request
from filelock import FileLock
from pandas_schema import Column, Schema
from pandas_schema.validation import MatchesPatternValidation, IsDistinctValidation

import snakemake
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import validate, min_version

import seq2science
from seq2science.util import (
    UniqueKeyLoader,
    parse_samples,
    samples2metadata,
    prep_filelock,
    url_is_alive,
    PickleDict,
    is_local,
    get_contrasts,
)


# get the cache and config dirs
CACHE_DIR = os.path.join(xdg.XDG_CACHE_HOME, "seq2science", seq2science.__version__)

# config.yaml(s)

# convert all config keys to lower case (upper case = user typo)
config = {k.lower(): v for k, v in config.items()}

# check config file for correct directory names
for key, value in config.items():
    if "_dir" in key:
        assert value not in [None, "", " "], f"\n{key} cannot be empty. For current directory, set '{key}: .'\n"
        # allow tilde as the first character, \w = all letters and numbers
        # ignore on Jenkins. it puts @ in the path
        assert (re.match("^[~\w_./-]*$", value[0]) and re.match("^[\w_./-]*$", value[1:])) or os.getcwd().startswith(
            "/var/lib/jenkins"
        ), (
            f"\nIn the config.yaml you set '{key}' to '{value}'. "
            + "Please use file paths that only contain letters, "
            + "numbers and any of the following symbols: underscores (_), periods (.), "
            + "and minuses (-). The first character may also be a tilde (~).\n"
        )
        config[key] = os.path.expanduser(value)

# make sure that difficult data-types (yaml objects) are in correct data format
for kw in ["aligner", "quantifier", "tpm2counts", "bam_sorter", "trimmer"]:
    if isinstance(config.get(kw, None), str):
        config[kw] = {config[kw]: {}}

# validate and complement the config dict
for schema in config_schemas:
    validate(config, schema=f"{config['rule_dir']}/../schemas/config/{schema}.schema.yaml")

# check if paired-end filename suffixes are lexicographically ordered
config["fqext"] = [config["fqext1"], config["fqext2"]]
assert sorted(config["fqext"])[0] == config["fqext1"], (
    "\nThe paired-end filename suffixes must be lexicographically ordered!\n"
    + f"Example suffixes: fqext1: R1, fqext2: R2\n"
    + f"Your suffixes:    fqext1: {config['fqext1']}, fqext2: {config['fqext2']}\n"
)

# read the config.yaml (not the profile)
with open(workflow.overwrite_configfiles[0], "r") as stream:
    # safe_load() with exception on duplicate keys
    user_config = yaml.load(stream, Loader=UniqueKeyLoader)

# make absolute paths, nest default dirs in result_dir and cut off trailing slashes
config["result_dir"] = re.split("\/$", os.path.abspath(config["result_dir"]))[0]

if not url_is_alive(config["samples"]):
    config["samples"] = os.path.abspath(config["samples"])

for key, value in config.items():
    if key.endswith("_dir"):
        if key in ["result_dir", "genome_dir", "rule_dir"] or key in user_config:
            value = os.path.abspath(value)
        else:
            value = os.path.abspath(os.path.join(config["result_dir"], value))
        config[key] = re.split("\/$", value)[0]


# samples.tsv


# read the samples.tsv file as all text, drop comment lines
try:
    samples = pd.read_csv(config["samples"], sep="\t", dtype="str", comment="#")
except Exception as e:
    logger.error("An error occurred while reading the samples.tsv:")
    column_error = re.search("Expected \d+ fields in line \d+, saw \d+", str(e.args))
    if column_error:
        digits = re.findall("\d+", column_error.group(0))
        logger.error(f"We found {digits[2]} columns on line {digits[1]}, but your header has only {digits[0]} columns...")
        show_header = True
        with open(config["samples"]) as s:
            for n, line in enumerate(s):
                if line.startswith("#"):
                    continue
                line = line.strip("\n").split("\t")
                if show_header:
                    logger.error(f"  header columns: {line}")
                    show_header = False
                    continue
                if n == int(digits[1]) - 1:
                    logger.error(f"  line {digits[1]} columns: {line}")
                    break
        logger.info("\n(if this was intentional, you can give this column an arbitrary name such as 'notes')")
        sys.exit(1)
    else:
        logger.error("")
        sys.exit(1)
samples.columns = samples.columns.str.strip()

assert all([col[0:7] not in ["Unnamed", ""] for col in samples]), (
    f"\nEncountered unnamed column in {config['samples']}.\n" + f"Column names: {str(', '.join(samples.columns))}.\n"
)

# use pandasschema for checking if samples file is filed out correctly
allowed_pattern = r"[A-Za-z0-9_.\-%]+"
distinct_columns = ["sample"]
if "descriptive_name" in samples.columns:
    distinct_columns.append("descriptive_name")

distinct_schema = Schema(
    [
        Column(
            col,
            [MatchesPatternValidation(allowed_pattern), IsDistinctValidation()]
            if col in distinct_columns
            else [MatchesPatternValidation(allowed_pattern)],
            allow_empty=True,
        )
        for col in samples.columns
    ]
)

errors = distinct_schema.validate(samples)

if len(errors):
    logger.error("\nThere are some issues with parsing the samples file:")
    for error in errors:
        logger.error(error)
    logger.error("")  # empty line
    sys.exit(1)

samples = parse_samples(samples, config)
if "technical_replicates" in samples:
    # check if replicate names are unique between assemblies
    r = samples[["assembly", "technical_replicates"]].drop_duplicates().set_index("technical_replicates")
    for replicate in r.index:
        assert len(r[r.index == replicate]) == 1, (
            "\nReplicate names must be different between assemblies.\n"
            + f"Replicate name '{replicate}' was found in assemblies {r[r.index == replicate]['assembly'].tolist()}."
        )

# check if sample, replicate and condition names are unique between the columns
for idx in samples.index:
    if "biological_replicates" in samples:
        assert idx not in samples["biological_replicates"].values, (
            f"sample names, conditions, and replicates can not overlap. "
            f"Sample {idx} can not also occur as a condition"
        )
    if "technical_replicates" in samples:
        assert idx not in samples["technical_replicates"].values, (
            f"sample names, conditions, and replicates can not overlap. "
            f"Sample {idx} can not also occur as a replicate"
        )

if "biological_replicates" in samples and "technical_replicates" in samples:
    for cond in samples["biological_replicates"]:
        assert cond not in samples["technical_replicates"].values, (
            f"sample names, conditions, and replicates can not overlap. "
            f"Condition {cond} can not also occur as a replicate"
        )

# validate samples file
for schema in sample_schemas:
    validate(samples, schema=f"{config['rule_dir']}/../schemas/samples/{schema}.schema.yaml")

sanitized_samples = copy.copy(samples)

samples = samples.set_index("sample")
samples.index = samples.index.map(str)


# check availability of assembly genomes and annotations


def get_workflow():
    # snakemake 6+ uses main_snakefile
    # return workflow.main_snakefile.split("/")[-2]
    return workflow.snakefile.split("/")[-2]


sequencing_protocol = (
    get_workflow()
    .replace("alignment", "Alignment")
    .replace("atac_seq", "ATAC-seq")
    .replace("chip_seq", "ChIP-seq")
    .replace("rna_seq", "RNA-seq")
    .replace("scatac_seq", "scATAC-seq")
)


if "assembly" in samples:
    # list assemblies that are used in this workflow
    used_assemblies = list(set(samples["assembly"]))

    # dictionary with which providers to use per genome
    providers = PickleDict(os.path.join(CACHE_DIR, "providers.p"))

    # split into local and remote assemblies, depending on if all required files can be found
    annotation_required = "rna_seq" in get_workflow() or config["aligner"] == "star"
    ftype = "annotation" if annotation_required else "genome"
    local_assemblies = [a for a in used_assemblies if is_local(a, ftype, config)]
    remote_assemblies = [a for a in used_assemblies if a not in local_assemblies]
    search_assemblies = [assembly for assembly in remote_assemblies if providers.get(assembly, {}).get(ftype) is None]

    # custom assemblies without provider (for local/remote annotations)
    if "scrna_seq" in get_workflow() and (
        "kite" in config["quantifier"].get("kallistobus", {}).get("ref", "")
        or config["quantifier"].get("citeseqcount", False)
    ):
        remote_assemblies = []
        search_assemblies = []
        
    # Check scRNA post-processing options
    if "scrna_seq" in get_workflow() and config.get("sc_preprocess",{}):
        if config["sc_preprocess"].get('run_sctk_qc',{}) and config["sc_preprocess"].get('export_sce_objects',{}):
            logger.error("Only one option is valid. Either select run_sctk_qc or export_sce_objects")
            sys.exit(1)
        
    # check if genomepy can download the required files
    if len(search_assemblies) > 0:
        providers.search(search_assemblies)

    # check if all remote assemblies can be downloaded
    if len(remote_assemblies) > 0:
        verbose = not config.get("no_config_log")
        providers.check(remote_assemblies, annotation_required, verbose)

    # trackhub

    # determine which assemblies (will) have an annotation
    _has_annot = dict()
    for assembly in local_assemblies:
        _has_annot[assembly] = is_local(assembly, "annotation", config)
        _has_annot[assembly + config.get("custom_assembly_suffix", "")] = _has_annot[assembly]
    for assembly in remote_assemblies:
        _has_annot[assembly] = bool(providers[assembly]["annotation"])
        _has_annot[assembly + config.get("custom_assembly_suffix", "")] = _has_annot[assembly]

    @lru_cache(maxsize=None)
    def has_annotation(assembly):
        """
        Returns True/False on whether or not the assembly has an annotation.
        """
        return _has_annot[assembly]

    # custom assemblies

    # control whether to custom extended assemblies
    if isinstance(config.get("custom_genome_extension"), str):
        config["custom_genome_extension"] = [config["custom_genome_extension"]]
    if isinstance(config.get("custom_annotation_extension"), str):
        config["custom_annotation_extension"] = [config["custom_annotation_extension"]]

    # custom assembly suffices
    modified = config.get("custom_genome_extension") or config.get("custom_annotation_extension")
    suffix = config["custom_assembly_suffix"] if modified else ""
    all_assemblies = [a + suffix for a in used_assemblies]

    def ori_assembly(assembly):
        """
        remove the extension suffix from an assembly if is was added.
        """
        return (
            assembly[: -len(config["custom_assembly_suffix"])]
            if assembly.endswith(config["custom_assembly_suffix"]) and modified
            else assembly
        )

    def custom_assembly(assembly):
        """
        add extension suffix to an assembly if is wasn't yet added.
        """
        return (
            assembly
            if assembly.endswith(config["custom_assembly_suffix"]) or not modified
            else (assembly + config["custom_assembly_suffix"])
        )

    # list of DESeq2 output files
    DE_contrasts = get_contrasts(samples, config, all_assemblies)

else:
    modified = False


# sample layouts


# check if a sample is single-end or paired end, and store it
if not config.get("no_config_log"):
    logger.info("Checking if samples are available online... This can take some time.")

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

            pickle.dump(
                {k: v for k, v in sampledict.items() if k.startswith(("ERR", "ERX", "SRR", "SRX", "GSM", "DRX", "DRR"))},
                open(pysradb_cache, "wb"),
            )

            # only keep samples for this run
            sampledict = {sample: values for sample, values in sampledict.items() if sample in all_samples}
        break
    except FileNotFoundError:
        time.sleep(1)
else:
    logger.error("There were some problems with locking the seq2science cache. Please try again in a bit.")
    sys.exit(1)

if not config.get("no_config_log"):
    logger.info("Done!\n\n")

# now check where to download which sample
ena_single_end = [
    run
    for values in sampledict.values()
    if (values["layout"] == "SINGLE") and values.get("ena_fastq_ftp") is not None
    for run in values["runs"]
]
ena_paired_end = [
    run
    for values in sampledict.values()
    if (values["layout"] == "PAIRED") and values.get("ena_fastq_ftp") is not None
    for run in values["runs"]
]
sra_single_end = [
    run
    for values in sampledict.values()
    if (values["layout"] == "SINGLE")
    for run in values.get("runs", [])
    if run not in ena_single_end
]
sra_paired_end = [
    run
    for values in sampledict.values()
    if (values["layout"] == "PAIRED")
    for run in values.get("runs", [])
    if run not in ena_paired_end
]

# get download link per run
run2download = dict()
for sample, values in sampledict.items():
    for run in values.get("runs", []):
        if values["ena_fastq_ftp"] and values["ena_fastq_ftp"][run]:
            if not (config.get("ascp_path") and config.get("ascp_key")):
                run2download[run] = [url.replace("era-fasp@fasp", "ftp") for url in values["ena_fastq_ftp"][run]]
            else:
                run2download[run] = values["ena_fastq_ftp"][run]


# TODO EXPLAIN!
rows = list()
merged_treps = set()
for sample, row in samples.iterrows():
    if len(sampledict[sample]["runs"]) > 0:
        for i, run in enumerate(sampledict[sample]["runs"]):
            row_dict = dict()
            rows.append(row_dict)

            row_dict["sample"] = run
            row_dict["technical_replicates"] = sample
            if "descriptive_name" in row:
                row_dict["descriptive_name"] = row.descriptive_name + f"_{i}"
            for col, value in [(col, value) for col, value in row.iteritems() if col not in ["sample", "descriptive_name", "technical_replicates"]]:
                row_dict[col] = value

            if i >= 1:
                merged_treps.add(sample)

    else:
        assert False

not_merged_treps = set(samples.index) - merged_treps
samples_extended_with_runs = pd.DataFrame(rows)
samples_extended_with_runs = samples_extended_with_runs.set_index("sample")
samples_extended_with_runs



# if samples are merged add the layout of the technical replicate to the config
failed_samples = dict()
if "technical_replicates" in samples_extended_with_runs:
    for replicate in samples_extended_with_runs.index:
        # replicate = samples.loc[sample, "technical_replicates"]
        if replicate not in sampledict:
            sampledict[replicate] = {"layout": sampledict[sample]["layout"]}
        elif sampledict[replicate]["layout"] != sampledict[sample]["layout"]:
            assembly = samples.loc[sample, "assembly"]
            treps = samples[(samples["assembly"] == assembly) & (samples["technical_replicates"] == replicate)].index
            failed_samples.setdefault(replicate, set()).update({trep for trep in treps})

if len(failed_samples):
    logger.error("Your technical replicates consist of a mix of single-end and paired-end samples!")
    logger.error("This is not supported.\n")

    for replicate, samples in failed_samples.items():
        logger.error(f"{replicate}:")
        for sample in samples:
            logger.error(f"\t{sample}: {sampledict[sample]['layout']}")
        logger.error("\n")
    sys.exit(1)


# workflow


def any_given(*args, prefix="", suffix=""):
    """
    returns a regex compatible string of all elements in the samples.tsv column given by the input
    """
    elements = []
    for column_name in args:
        if column_name == "sample":
            elements.extend(samples.index)
        elif column_name in samples:
            elements.extend(samples[column_name])

    elements = [prefix + element + suffix for element in elements if isinstance(element, str)]
    return "|".join(set(elements))


# set global wildcard constraints (see workflow._wildcard_constraints)
sample_constraints = ["sample"]


wildcard_constraints:
    sorting="coordinate|queryname",
    sorter="samtools|sambamba",


if "assembly" in samples:

    wildcard_constraints:
        raw_assembly=any_given("assembly"),
        assembly=any_given("assembly", suffix=config["custom_assembly_suffix"] if modified else ""),


if "technical_replicates" in samples:
    sample_constraints.append("technical_replicates")

    wildcard_constraints:
        replicate=any_given("technical_replicates"),


if "biological_replicates" in samples:
    sample_constraints.append("biological_replicates")

    wildcard_constraints:
        condition=any_given("biological_replicates"),


if "control" in samples:
    sample_constraints.append("control")


wildcard_constraints:
    sample=any_given(*sample_constraints),


# make sure the snakemake version corresponds to version in environment
if snakemake.__version__ != "5.18.1":
    raise WorkflowError(
        f"Expecting Snakemake version 5.18.1 (you are currently using {snakemake.__version__}"
    )

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
                    try:
                        response = requests.get(f"https://genome.ucsc.edu/cgi-bin/hgGateway", allow_redirects=True)
                    except:
                        logger.error("There seems to be some problems with connecting to UCSC, try again in some time")
                        assert False
                    if not response.ok:
                        logger.error("Make sure you are connected to the internet")
                        assert False

                    with urllib.request.urlopen("https://api.genome.ucsc.edu/list/ucscGenomes") as url:
                        data = json.loads(url.read().decode())["ucscGenomes"]

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
        sys.exit(1)


def config_rerun_parser(configdict):
    """
    Makes a copy of the configdict with certain keys removed since they are not relevant
    to rerunning inference.
    """
    configdict = copy.copy(configdict)
    configdict.pop("no_config_log", None)
    configdict.pop("cli_call", None)
    configdict.pop("cores", None)

    return configdict
