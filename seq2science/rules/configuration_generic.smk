"""
all logic not related to any specific workflows should be here.
"""
container: "docker://quay.io/biocontainers/seq2science:0.7.2--pypyhdfd78af_0"

import os.path
import pickle
import re
import time
import copy
import json
import requests
import yaml

import xdg
import pandas as pd
import urllib.request
from filelock import FileLock
from pandas_schema import Column, Schema
from pandas_schema.validation import MatchesPatternValidation, IsDistinctValidation

from snakemake.dag import DAG
from snakemake.logging import logger
from snakemake.utils import validate

import seq2science
from seq2science.util import (
    UniqueKeyLoader,
    dense_samples,
    samples2metadata,
    prep_filelock,
    url_is_alive,
    PickleDict,
    is_local,
    get_contrasts,
    assert_versions,
)


# monkeypatch snakemake to not see which files are changed
# since we fix that ourselves
DAG.warn_about_changes = lambda *args: None

# get the cache and config dirs
CACHE_DIR = os.path.join(xdg.XDG_CACHE_HOME, "seq2science", seq2science.__version__)


def parse_config(cfg):
    # convert all config keys to lower case (upper case = user typo)
    cfg = {k.lower(): v for k, v in cfg.items()}
    
    # check config file for correct directory names
    for key, value in cfg.items():
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
            cfg[key] = os.path.expanduser(value)
    
    # make sure that difficult data-types (yaml objects) are in correct data format
    for kw in ["aligner", "quantifier", "tpm2counts", "bam_sorter", "trimmer"]:
        if isinstance(cfg.get(kw, None), str):
            cfg[kw] = {cfg[kw]: {}}
    
    # validate and complement the config dict
    for schema in CONFIG_SCHEMAS:
        validate(cfg, schema=f"{cfg['rule_dir']}/../schemas/config/{schema}.schema.yaml")
    
    # check if paired-end filename suffixes are lexicographically ordered
    cfg["fqext"] = [cfg["fqext1"], cfg["fqext2"]]
    assert sorted(cfg["fqext"])[0] == cfg["fqext1"], (
        "\nThe paired-end filename suffixes must be lexicographically ordered!\n"
        + f"Example suffixes: fqext1: R1, fqext2: R2\n"
        + f"Your suffixes:    fqext1: {cfg['fqext1']}, fqext2: {cfg['fqext2']}\n"
    )
    
    # read the config.yaml (not the profile)
    with open(workflow.overwrite_configfiles[0], "r") as stream:
        # safe_load() with exception on duplicate keys
        user_config = yaml.load(stream, Loader=UniqueKeyLoader)
    
    # make absolute paths, nest default dirs in result_dir and cut off trailing slashes
    cfg["result_dir"] = re.split("/$", os.path.abspath(cfg["result_dir"]))[0]
    
    if not url_is_alive(cfg["samples"]):
        cfg["samples"] = os.path.abspath(cfg["samples"])
    
    for key, value in cfg.items():
        if key.endswith("_dir"):
            if key in ["result_dir", "genome_dir", "rule_dir"] or key in user_config:
                value = os.path.abspath(value)
            else:
                value = os.path.abspath(os.path.join(cfg["result_dir"], value))
            cfg[key] = re.split("/$", value)[0]
    
    # check that the versions in the requirements.yaml match the installed versions
    assert_versions(cfg["rule_dir"])
    
    return cfg

config = parse_config(config)  # overwrite the existing global


# samples.tsv


def parse_samples():
    # read the samples.tsv file as all text, drop comment lines
    samples_df = pd.read_csv(
        config["samples"], sep="\t", dtype="str",
        comment="#", index_col=False, header=0,
    )

    def parsing_error(messages: str or list):
        if isinstance(messages, str):
            messages = [messages]
        logger.error("\nThere are some issues with the samples file:")
        for message in messages:
            logger.error(message)
        logger.error("")  # empty line
        os._exit(1)  # noqa

    # check column names
    samples_df.columns = samples_df.columns.str.strip()
    errors = []
    for col in samples_df.columns:
        violations = re.findall(r"[^A-Za-z0-9_.\-%]", col)
        if len(violations):
            n = len(violations)
            violations = '", "'.join(violations)
            errors.append(
                f'The column "{col}" contains {n} forbidden symbol{"" if n == 1 else "s"}: "{violations}"'
            )
    if len(errors):
        parsing_error(errors)
    if any([col[0:7] in ["Unnamed", ""] for col in samples_df]):
        columns = '", "'.join(samples_df.columns)
        parsing_error(f'Encountered unnamed column(s): "{columns}"')

    # check dataframe shape
    # (rows longer than the number of columns are truncated by pandas)
    with open(config["samples"]) as s:
        for n, line in enumerate(s):
            if line.startswith("#"):
                continue
            line = line.split("\t")
            if len(line) > len(samples_df.columns):
                errors.append(
                    f"Line {n} contains {len(line)} fields, but there are only {len(samples_df.columns)} column names!"
                )
    if len(errors):
        cols = '", "'.join(samples_df.columns)
        errors.append(f'Columns names: "{cols}"')
        errors.append('(if this was intentional, you can give columns arbitrary names such as "notes")')
        parsing_error(errors)

    # use pandasschema for checking if samples file is filed out correctly
    allowed_pattern = r"^[A-Za-z0-9_.\-%]+$"
    spellcheck_columns = ["sample"]
    for col in ["assembly", "technical_replicates", "descriptive_name", "biological_replicates", "control"]:
        if col in samples_df.columns:
            spellcheck_columns.append(col)
    schema = Schema(
        [Column(col, [MatchesPatternValidation(allowed_pattern)], allow_empty=True) for col in spellcheck_columns]
    )
    errors = schema.validate(samples_df, spellcheck_columns)

    distinct_columns = ["sample"]
    if "descriptive_name" in samples_df.columns:
        distinct_columns.append("descriptive_name")
    schema = Schema(
        [Column(col, [IsDistinctValidation()], allow_empty=True) for col in distinct_columns]
    )
    errors += schema.validate(samples_df, distinct_columns)
    
    if len(errors):
        logger.error("\nThere are some issues with the samples file:")
        spellcheck_error = f'does not match the pattern "{allowed_pattern}"'
        for error in errors:
            error = str(error)
            if error.endswith(spellcheck_error):
                violations = re.findall(r"[^A-Za-z0-9_.\-%]", error.split('"')[3])
                n = len(violations)
                violations = '", "'.join(violations)
                error = error.replace(
                    spellcheck_error,
                    f'contain {n} forbidden symbol{"" if n==1 else "s"}: "{violations}"'
                )
            logger.error(error)
        logger.error("")  # empty line
        os._exit(1)  # noqa

    # remove unused columns, and populate empty cells in used columns
    samples_df = dense_samples(samples_df, config)

    # check if replicate names are unique between assemblies
    if "technical_replicates" in samples_df:
        r = samples_df[["assembly", "technical_replicates"]].drop_duplicates().set_index("technical_replicates")
        for replicate in r.index:
            assert len(r[r.index == replicate]) == 1, (
                "\nReplicate names must be different between assemblies.\n"
                + f"Replicate name '{replicate}' was found in assemblies {r[r.index == replicate]['assembly'].tolist()}."
            )
    
    # check if sample, replicate and condition names are unique between the columns
    for idx in samples_df.index:
        if "biological_replicates" in samples_df:
            assert idx not in samples_df["biological_replicates"].values, (
                f"sample names, conditions, and replicates can not overlap. "
                f"Sample {idx} can not also occur as a condition"
            )
        if "technical_replicates" in samples_df:
            assert idx not in samples_df["technical_replicates"].values, (
                f"sample names, conditions, and replicates can not overlap. "
                f"Sample {idx} can not also occur as a replicate"
            )
    
    if "biological_replicates" in samples_df and "technical_replicates" in samples_df:
        for cond in samples_df["biological_replicates"]:
            assert cond not in samples_df["technical_replicates"].values, (
                f"sample names, conditions, and replicates can not overlap. "
                f"Condition {cond} can not also occur as a replicate"
            )
    
    # validate samples file
    for schema in SAMPLE_SCHEMAS:
        validate(samples_df, schema=f"{config['rule_dir']}/../schemas/samples/{schema}.schema.yaml")
    
    sanitized_samples_df = copy.copy(samples_df)
    
    samples_df = samples_df.set_index("sample")
    samples_df.index = samples_df.index.map(str)
    
    return samples_df, sanitized_samples_df

samples, sanitized_samples = parse_samples()

# all samples (including controls)
ALL_SAMPLES = [sample for sample in samples.index]
if "control" in samples.columns:
    ALL_SAMPLES.extend(set(samples["control"].dropna()))


# technical and biological replicates
if "assembly" in samples:
    def get_replicates():
        # dataframe with all technical replicates collapsed
        cols = ["sample", "assembly"]
        subset = ["sample", "assembly"]
        if "technical_replicates" in samples:
            cols = ["technical_replicates", "assembly"]
            subset = ["technical_replicates", "assembly"]
        if "biological_replicates" in samples:
            cols.append("biological_replicates")
            subset.append("biological_replicates")
        if "control" in samples:
            cols.append("control")
        if "colors" in samples:
            cols.append("colors")
        treps = samples.reset_index()[cols].drop_duplicates(subset=subset).set_index(cols[0])
        assert treps.index.is_unique, f"duplicate value found in treps: {treps[treps.index.duplicated()]}"

        # dataframe with all replicates collapsed
        breps = treps
        if "biological_replicates" in treps:
            breps = treps.reset_index(drop=True).drop_duplicates(subset=subset[1:]).set_index("biological_replicates")

        # treps that came from a merge
        merged_treps = [trep for trep in treps.index if trep not in samples.index]

        # make a dict that returns the treps that belong to a brep
        treps_from_brep = dict()
        if "biological_replicates" in treps:
            for brep, row in breps.iterrows():
                assembly = row["assembly"]
                treps_from_brep[(brep, assembly)] = list(
                    treps[(treps["assembly"] == assembly) & (treps["biological_replicates"] == brep)].index
                )
                treps_from_brep[(brep, assembly + config.get("custom_assembly_suffix",""))] = list(
                    treps[(treps["assembly"] == assembly) & (treps["biological_replicates"] == brep)].index
                )
        else:
            for brep, row in breps.iterrows():
                assembly = row["assembly"]
                treps_from_brep[(brep, assembly)] = [brep]
                treps_from_brep[(brep, assembly + config.get("custom_assembly_suffix",""))] = [brep]

        # and vice versa
        brep_from_treps = dict()
        for (brep, _assembly), _treps in treps_from_brep.items():
            brep_from_treps.update({trep: brep for trep in _treps})

        return treps, merged_treps, breps, treps_from_brep, brep_from_treps

    treps, MERGED_TREPS, breps, TREPS_FROM_BREP, BREP_FROM_TREPS = get_replicates()


def rep_to_descriptive(rep, brep=False):
    """
    Return the descriptive name for a replicate.
    """
    if brep and "biological_replicates" in samples:
        rep = samples[samples.biological_replicates == rep].biological_replicates[0]
    elif "descriptive_name" in samples:
        col = samples.index
        if "technical_replicates" in samples:
            col = samples.technical_replicates
        rep = samples[col == rep].descriptive_name[0]
    return rep


# check availability of assembly genomes and annotations

WORKFLOW = workflow.main_snakefile.split("/")[-2]

SEQUENCING_PROTOCOL = (
    WORKFLOW
    .replace("alignment", "Alignment")
    .replace("atac_seq", "ATAC-seq")
    .replace("chip_seq", "ChIP-seq")
    .replace("rna_seq", "RNA-seq")
    .replace("scatac_seq", "scATAC-seq")
)

modified = False
if "assembly" in samples:
    def parse_assemblies():
        # list assemblies that are used in this workflow
        used_assemblies = list(set(samples["assembly"]))
        if "motif2factors_reference" in config and config["run_gimme_maelstrom"]:
            _used_assemblies = used_assemblies + config["motif2factors_reference"]
    
        # dictionary with which providers to use per genome
        providers = PickleDict(os.path.join(CACHE_DIR, "providers.p"))
    
        # split into local and remote assemblies, depending on if all required files can be found
        annotation_required = "rna_seq" in WORKFLOW or config["aligner"] == "star"
        ftype = "annotation" if annotation_required else "genome"
        local_assemblies = [a for a in _used_assemblies if is_local(a, ftype, config)]
        remote_assemblies = [a for a in _used_assemblies if a not in local_assemblies]
        search_assemblies = [assembly for assembly in remote_assemblies if providers.get(assembly, {}).get(ftype) is None]
    
        # custom assemblies without provider (for local/remote annotations)
        if "scrna_seq" in WORKFLOW and (
            "kite" in config["quantifier"].get("kallistobus", {}).get("ref", "")
            or config["quantifier"].get("citeseqcount", False)
        ):
            remote_assemblies = []
            search_assemblies = []
            
        # Check scRNA post-processing options
        if "scrna_seq" in WORKFLOW and config.get("sc_preprocess",{}):
            if config["sc_preprocess"].get('run_sctk_qc',{}) and config["sc_preprocess"].get('export_sce_objects',{}):
                logger.error("Only one option is valid. Either select run_sctk_qc or export_sce_objects")
                os._exit(1)  # noqa
            
        # check if genomepy can download the required files
        if len(search_assemblies) > 0:
            providers.search(search_assemblies)
    
        # check if all remote assemblies can be downloaded
        if len(remote_assemblies) > 0:
            verbose = not config.get("no_config_log")
            providers.check(remote_assemblies, annotation_required, verbose)
    
        # trackhub
    
        # determine which assemblies (will) have an annotation
        has_annotation = dict()
        for assembly in local_assemblies:
            has_annotation[assembly] = is_local(assembly,"annotation", config)
            if config.get("custom_assembly_suffix", False):
                has_annotation[assembly + config["custom_assembly_suffix"]] = has_annotation[assembly]
        for assembly in remote_assemblies:
            has_annotation[assembly] = bool(providers[assembly]["annotation"])
            if config.get("custom_assembly_suffix", False):
                has_annotation[assembly + config["custom_assembly_suffix"]] = has_annotation[assembly]
    
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
            remove the extension suffix from an assembly if it was added.
            """
            if not modified or not assembly.endswith(config["custom_assembly_suffix"]):
                return assembly
            return assembly[: -len(config["custom_assembly_suffix"])]

        ori_assemblies    = {a: ori_assembly(a) for a in all_assemblies}
        custom_assemblies = {ori_assembly(a): a for a in all_assemblies}
    
        # list of DESeq2 output files
        de_contrasts = get_contrasts(samples, config, all_assemblies)
        
        return providers, has_annotation, all_assemblies, modified, ori_assemblies, custom_assemblies, de_contrasts

    PROVIDERS, HAS_ANNOTATION, ALL_ASSEMBLIES, modified, ORI_ASSEMBLIES, CUSTOM_ASSEMBLIES, DE_CONTRASTS = parse_assemblies()

# sample layouts

def parse_pysradb():
    # check if a sample is single-end or paired end, and store it
    if not config.get("no_config_log"):
        logger.info("Checking if samples are available online... This can take some time.")

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

                missing_samples = [sample for sample in ALL_SAMPLES if sample not in sampledict.keys()]
                if len(missing_samples) > 0:
                    sampledict.update(samples2metadata(missing_samples, config, logger))

                pickle.dump(
                    {k: v for k, v in sampledict.items() if k.startswith(("ERR", "ERX", "SRR", "SRX", "GSM", "DRX", "DRR"))},
                    open(pysradb_cache, "wb"),
                )

                # only keep samples for this run
                sampledict = {sample: values for sample, values in sampledict.items() if sample in ALL_SAMPLES}
            break
        except FileNotFoundError:
            time.sleep(1)
    else:
        logger.error("There were some problems with locking the seq2science cache. Please try again in a bit.")
        os._exit(1)  # noqa

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

    # if samples are merged add the layout of the technical replicate to the config
    failed_samples = dict()
    if "technical_replicates" in samples:
        for sample in samples.index:
            replicate = samples.loc[sample, "technical_replicates"]
            if replicate not in sampledict:
                sampledict[replicate] = {"layout": sampledict[sample]["layout"]}
            elif sampledict[replicate]["layout"] != sampledict[sample]["layout"]:
                assembly = samples.loc[sample, "assembly"]
                treps = samples[(samples["assembly"] == assembly) & (samples["technical_replicates"] == replicate)].index
                failed_samples.setdefault(replicate, set()).update({trep for trep in treps})

    if len(failed_samples):
        logger.error("Your technical replicates consist of a mix of single-end and paired-end samples!")
        logger.error("This is not supported.\n")

        for replicate, failed_samples in failed_samples.items():
            logger.error(f"{replicate}:")
            for sample in failed_samples:
                logger.error(f"\t{sample}: {sampledict[sample]['layout']}")
            logger.error("\n")
        os._exit(1)  # noqa

    return sampledict, ena_single_end, ena_paired_end, sra_single_end, sra_paired_end, run2download, pysradb_cache_lock

SAMPLEDICT, ENA_SINGLE_END, ENA_PAIRED_END, SRA_SINGLE_END, SRA_PAIRED_END, RUN2DOWNLOAD, PYSRADB_CACHE_LOCK = parse_pysradb()


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
            if (column_name == "assembly") and ("motif2factors_reference" in config):
                elements.extend(config["motif2factors_reference"])
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


if config.get("create_trackhub"):
    def parse_ucsc():
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
            os._exit(1)  # noqa

        return ucsc_assemblies

    # record which assembly trackhubs are found on UCSC
    UCSC_ASSEMBLIES = parse_ucsc()


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


if 'alignment_general' in CONFIG_SCHEMAS:
    # Depending on the config settings, one of three rules in bam_cleaning.smk
    # will output the final bams: mark_duplicates, sieve_bam or samtools_sort
    FINAL_BAM = f"{config['final_bam_dir']}/{{assembly}}-{{sample}}.samtools-coordinate.bam"
    FINAL_BAI = f"{FINAL_BAM}.bai"
