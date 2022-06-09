"""
Utility functions for seq2science
"""
import contextlib
import re
import os
import sys
import glob
import time
import colorsys
import pickle
import urllib.request
import yaml
from io import StringIO
from typing import List
from math import ceil, floor

from filelock import FileLock
import genomepy
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pysradb
from snakemake.logging import logger

# default colors in matplotlib. Order dictates the priority.
DEFAULT_COLOR_DICTS = [mcolors.BASE_COLORS, mcolors.TABLEAU_COLORS, mcolors.CSS4_COLORS, mcolors.XKCD_COLORS]


class UniqueKeyLoader(yaml.SafeLoader):
    """
    YAML loader with duplicate key checking
    source: https://stackoverflow.com/questions/33490870/parsing-yaml-in-python-detect-duplicated-keys
    """

    def construct_mapping(self, node, deep=False):
        mapping = []
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep).lower()
            if key in mapping:
                logger.error(f"Duplicate key found in the config.yaml: {key}\n")
                os._exit(1)  # noqa
            mapping.append(key)
        return super().construct_mapping(node, deep)


def _sample_to_idxs(df: pd.DataFrame, sample: str) -> List[int]:
    """
    Get a list of index/indices that belong to a run in a pysradb dataframe
    """
    if sample.startswith(("SRR", "DRR", "ERR")):
        idxs = df.index[df.run_accession == sample].tolist()
        assert len(idxs) == 1, f"sample {sample} with idxs: {idxs}"
    elif sample.startswith(("SRX", "ERX", "DRX")):
        idxs = df.index[df.experiment_accession == sample].tolist()
        assert len(idxs) >= 1, len(idxs)
    else:
        assert False, (
            f"sample {sample} not a run, this should not be able to happen!" f" Please make an issue about this!"
        )
    return idxs


def dense_samples(samples: pd.DataFrame, config: dict) -> pd.DataFrame:
    """
    for each functional column, if found in samples.tsv:
    1) if it is incomplete, fill the blanks with replicate/sample names
       (sample names if replicates are not found/applicable)
    2) drop column if it provides no information
       (renamed in case it's used in a DE contrast)
    """
    if "technical_replicates" in samples:
        samples["technical_replicates"] = samples["technical_replicates"].mask(pd.isnull, samples["sample"])
        if (
                len(samples["technical_replicates"].unique()) == len(samples["sample"].unique())
                or config.get("technical_replicates") == "keep"
        ):
            samples.rename(columns={"technical_replicates": "_trep"}, inplace=True)
    col = "technical_replicates" if "technical_replicates" in samples else "sample"
    if "biological_replicates" in samples:
        samples["biological_replicates"] = samples["biological_replicates"].mask(pd.isnull, samples[col])
        if (
                len(samples["biological_replicates"].unique()) == len(samples[col].unique())
                or config.get("biological_replicates") == "keep"
        ):
            samples.rename(columns={"biological_replicates": "_brep"}, inplace=True)
    if "descriptive_name" in samples:
        samples["descriptive_name"] = samples["descriptive_name"].mask(pd.isnull, samples[col])
        if samples["descriptive_name"].to_list() == samples[col].to_list():
            samples.rename(columns={"descriptive_name": "_dname"}, inplace=True)
    if "strandedness" in samples:
        samples["strandedness"] = samples["strandedness"].mask(pd.isnull, "nan")
        if config.get("ignore_strandedness", True) or not any(
                [field in list(samples["strandedness"]) for field in ["yes", "forward", "reverse", "no"]]
        ):
            samples = samples.drop(columns=["strandedness"])
    if "colors" in samples:
        if config.get("create_trackhub", False):
            samples["colors"] = samples["colors"].mask(pd.isnull, "0,0,0")  # nan -> black
            samples["colors"] = [color_parser(c) for c in samples["colors"]]  # convert input to HSV color
        else:
            samples = samples.drop(columns=["colors"])
    return samples


def samples2metadata_local(samples: List[str], config: dict, logger) -> dict:
    """
    (try to) get the metadata of local samples
    """
    SAMPLEDICT = dict()
    for sample in samples:
        local_fastqs = glob.glob(os.path.join(config["fastq_dir"], f'{sample}*{config["fqsuffix"]}*.gz'))
        if len(local_fastqs) == 1:
            SAMPLEDICT[sample] = dict()
            SAMPLEDICT[sample]["layout"] = "SINGLE"
        elif (
            len(local_fastqs) == 2
            and any([config["fqext1"] in os.path.basename(f) for f in local_fastqs])
            and any([config["fqext2"] in os.path.basename(f) for f in local_fastqs])
        ):
            SAMPLEDICT[sample] = dict()
            SAMPLEDICT[sample]["layout"] = "PAIRED"
        elif sample.startswith(("GSM", "DRX", "ERX", "SRX", "DRR", "ERR", "SRR")):
            continue
        else:
            extend_msg = ""
            if len(local_fastqs) > 2:
                extend_msg = (
                    f"We found too many files matching ({len(local_fastqs)}) "
                    "and could not distinguish them:\n" + ", ".join([os.path.basename(f) for f in local_fastqs]) + ".\n"
                )

            logger.error(
                f"\nsample {sample} was not found..\n"
                f"We checked directory '{config['fastq_dir']}' "
                f"for gzipped files starting with '{sample}' and containing '{config['fqsuffix']}'.\n"
                + extend_msg
                + f"Since the sample did not start with either GSM, SRX, SRR, ERR, and DRR we "
                f"couldn't find it online..\n"
            )
            os._exit(1)  # noqa

    return SAMPLEDICT


def samples2metadata_sra(samples: List[str], logger) -> dict:
    """
    Get the required info to continue a seq2science run from a list of samples.

    - If a sample already exists locally, we only want to know if it is paired-end or single-end.
    - If a sample does not exist locally
      - find its corresponding SRX number and all runs that belong to it,
      - check if they all have the same layout, if not, crash
      - see if we can download the runs from ena

    output:
        dict(
            "GSM1234": {"layout": "PAIRED",
                         "runs": ["SRR1234", "SRR4321"],
                         "ena_fastq_ftp": {...},

            "SRR5678": {"layout": "SINGLE",
                        "runs": ["SRR5678"],
                        ena_fastq_ftp: None,
            ...
        )
    """
    # start with empty dictionary which we fill out later
    SAMPLEDICT = {sample: dict() for sample in samples}

    # only continue with public samples
    db_sra = pysradb.SRAweb()

    # all samples that are on GEO (GSM numbers), must first be converted to SRA ids (SRX numbers)
    geo_samples = [sample for sample in samples if sample.startswith("GSM")]

    # in sample2clean we store the (potential GEO) sample name in a SRA compliant name
    if len(geo_samples):
        try:
            df_geo = db_sra.gsm_to_srx(geo_samples)
        except:
            logger.error(
                "We had trouble querying the SRA. This probably means that the SRA was unresponsive, and their servers "
                "are overloaded or slow. Please try again in a bit...\n"
                "Another possible option is that you try to access samples that do not exist or are protected, and "
                "seq2science does not support downloading those..\n\n"
            )
            os._exit(1)  # noqa

        sample2clean = dict(zip(df_geo.experiment_alias, df_geo.experiment_accession))
    else:
        sample2clean = dict()

    # now add the already SRA compliant names with a reference to itself
    sample2clean.update({sample: sample for sample in samples if sample not in geo_samples})

    # check our samples on sra
    try:
        df_sra = db_sra.sra_metadata(list(sample2clean.values()), detailed=True)
    except:
        logger.error(
            "We had trouble querying the SRA. This probably means that the SRA was unresponsive, and their servers "
            "are overloaded or slow. Please try again in a bit...\n"
            "Another possible option is that you try to access samples that do not exist or are protected, and "
            "seq2science does not support downloading those..\n\n"
        )
        os._exit(1)  # noqa

    # keep track of not-supported samples
    not_supported_formats = ["ABI_SOLID"]
    not_supported_samples = []

    for sample, clean in sample2clean.items():
        # table indices
        idxs = _sample_to_idxs(df_sra, clean)

        # get all runs that belong to the sample
        runs = df_sra.loc[idxs].run_accession.tolist()
        assert len(runs) >= 1
        SAMPLEDICT[sample]["runs"] = runs

        # check if sample is from a supported format
        for bad_format in not_supported_formats:
            for real_format in df_sra.loc[idxs].instrument_model_desc.tolist():
                if real_format == bad_format:
                    not_supported_samples.append(sample)

        # get the layout
        layout = df_sra.loc[idxs].library_layout.tolist()
        assert len(set(layout)) == 1, f"sample {sample} consists of mixed layouts, bad!"
        assert layout[0] in ["PAIRED", "SINGLE"], f"sample {sample} is an unclear layout, bad!"
        SAMPLEDICT[sample]["layout"] = layout[0]

        # get the ena url
        SAMPLEDICT[sample]["ena_fastq_ftp"] = dict()
        for run in runs:
            if layout[0] == "SINGLE":
                SAMPLEDICT[sample]["ena_fastq_ftp"][run] = df_sra[df_sra.run_accession == run].ena_fastq_ftp.tolist()
            elif layout[0] == "PAIRED":
                SAMPLEDICT[sample]["ena_fastq_ftp"][run] = (
                    df_sra[df_sra.run_accession == run].ena_fastq_ftp_1.tolist()
                    + df_sra[df_sra.run_accession == run].ena_fastq_ftp_2.tolist()
                )

        # if any run from a sample is not found on ENA, better be safe, and assume that sample as a whole is not on ENA
        if any([any(pd.isna(urls)) for urls in SAMPLEDICT[sample]["ena_fastq_ftp"].values()]):
            SAMPLEDICT[sample]["ena_fastq_ftp"] = None

    # now report single message for all sample(s) that are from a sequencing platform that is not supported
    assert len(not_supported_samples) == 0, (
        f'Sample(s) {", ".join(not_supported_samples)} are not supported by seq2science. Samples that are one of '
        f'these formats; [{", ".join(not_supported_formats)}] are not supported.'
    )

    return SAMPLEDICT


def samples2metadata(samples: List[str], config: dict, logger) -> dict:
    local_samples = samples2metadata_local(samples, config, logger)
    public_samples = [sample for sample in samples if sample not in local_samples.keys()]

    if len(public_samples) == 0:
        return local_samples

    # chop public samples into smaller chunks, doing large queries results into
    # pysradb decode errors..
    chunksize = 100
    chunked_public = [public_samples[i : i + chunksize] for i in range(0, len(public_samples), chunksize)]
    sra_samples = dict()
    for chunk in chunked_public:
        sra_samples.update(samples2metadata_sra(chunk, logger))
        # just to be sure sleep in between to not go over our API limit
        time.sleep(1)

    return {**local_samples, **sra_samples}


def sieve_bam(configdict):
    """
    helper function to check whether to use rule sieve_bam
    """
    return (
        configdict.get("min_mapping_quality", 0) > 0
        or configdict.get("tn5_shift", False)
        or configdict.get("remove_blacklist", False)
        or configdict.get("filter_on_size", False)
        or configdict.get("remove_mito", False)
        or configdict.get("supsample", -1) > 0
    )


def parse_contrast(contrast, samples, check=True):
    """
    Extract contrast batch and column, target and reference groups from a DE contrast design.
    Check for contrast validity if check = True.

    If "all" is in the contrast groups, it is always assumed to be the target
    (and expanded to mean all groups in the column in function `get_contrasts()`).

    Accepts a string containing a DESeq2 contrast design.

    Returns
        batch: the batch column, or None
        column: the contrast column
        target: the group of interest
        reference: the control group
    """
    # clean contrast
    de_contrast = contrast.strip().replace(" ", "").replace("~", "")

    # parse batch effect
    batch = None
    if "+" in de_contrast:
        batch, de_contrast = de_contrast.split("+")
        if len(de_contrast) > 1:
            ValueError(f"DE contrast {contrast} can only contain a '+' to denote the batch effect column.")

    # parse groups
    target, reference = de_contrast.split("_")[-2:]

    # parse column
    n = de_contrast.find(f"_{target}_{reference}")
    column = de_contrast[:n]

    # "all" is never the reference
    if reference == "all":
        reference = target
        target = "all"

    if check:
        # check if columns exists and are valid
        valid_columns = [col for col in samples.columns if col not in ["sample", "assembly"]]
        # columns that may have been dropped, if so, these backups have been saved
        backup_columns = {"technical_replicates": "_trep", "biological_replicates": "_brep", "descriptive_name": "_dname"}
        for col in [batch, column]:
            if col:
                assert col in valid_columns + list(backup_columns.keys()), (
                    f'\nIn DESeq2 contrast design "{contrast}", '
                    f'column "{col}" does not match any valid column name.\n'
                )
        # check if group is element of column
        c = column if column in samples else backup_columns[column]  # column/backup column
        all_groups = set(samples[c].dropna().astype(str))
        for group in [target, reference]:
            if group != "all":
                assert group in all_groups, (
                    f'\nIn DESeq2 contrast design "{contrast}", ' f"group {group} cannot be found in column {column}.\n"
                )

    return batch, column, target, reference


def expand_contrasts(samples: pd.DataFrame, contrasts: list or str) -> list:
    """
    splits contrasts that contain multiple comparisons
    """
    if contrasts is None:
        return []
    if isinstance(contrasts, str):
        contrasts = [contrasts]

    new_contrasts = []
    for contrast in contrasts:
        batch, column, target, reference = parse_contrast(contrast, samples, check=False)

        if target == "all":
            # all vs 1 comparison ("all vs A")
            targets = set(samples[column].dropna().astype(str))
            targets.remove(reference)
        else:
            # 1 vs 1 comparison ("A vs B")
            targets = [target]

        for target in targets:
            new_contrast = f"{column}_{target}_{reference}"
            if batch:
                new_contrast = f"{batch}+{new_contrast}"
            new_contrasts.append(new_contrast)

    # get unique elements
    new_contrasts = list(set(new_contrasts))
    return new_contrasts


def get_contrasts(samples: pd.DataFrame, config: dict, ALL_ASSEMBLIES: list) -> list:
    """
    list all diffexp.tsv files we expect
    """
    contrasts = expand_contrasts(samples, config.get("contrasts"))
    all_contrasts = list()
    for de_contrast in contrasts:
        # parse groups
        target, reference = de_contrast.split("_")[-2:]
        column = de_contrast[: de_contrast.find(f"_{target}_{reference}")]

        # parse column
        if "+" in column:
            column = column.split("+")[1]

        if column not in samples:
            backup_columns = {
                "technical_replicates": "_trep",  # is trep technically possible? You need multiple reps right?
                "biological_replicates": "_brep",
                "descriptive_name": "_dname",
            }
            column = backup_columns[column]

        for assembly in ALL_ASSEMBLIES:
            groups = set(samples[samples.assembly == assembly][column].to_list())
            if target in groups and reference in groups:
                all_contrasts.append(f"{config['deseq2_dir']}/{assembly}-{de_contrast}.diffexp.tsv")
    return all_contrasts


def url_is_alive(url):
    """
    Checks that a given URL is reachable.
    https://gist.github.com/dehowell/884204
    """
    for i in range(3):
        try:
            request = urllib.request.Request(url)
            request.get_method = lambda: "HEAD"

            urllib.request.urlopen(request, timeout=5)
            return True
        except:
            continue
    return False


def get_bustools_rid(params):
    """
    Extract the position of the fastq containig reads from the bustools -x argument.
    The read_id is the first pos of the last triplet in the bc:umi:read string or hard-coded
    for short-hand syntax.
    In: -x 10xv3 -> read_id=1
    In: -x 0,0,16:0,16,26:1,0,0 -> read_id=1
    """
    kb_tech_dict = {
        "10xv2": 1,
        "10xv3": 1,
        "celseq": 1,
        "celseq2": 1,
        "dropseq": 1,
        "scrubseq": 1,
        "indropsv1": 1,
        "indropsv2": 0,
    }
    # Check for occurence of short-hand tech
    bus_regex = "(?<!\S)([0-1],\d*,\d*:){2}([0-1],0,0)(?!\S)"
    bus_regex_short = "(?i)\\b(10XV2|10XV3|CELSEQ|CELSEQ2|DROPSEQ|SCRUBSEQ|INDROPSV1|INDROPSV2)\\b"

    if re.search(bus_regex, params) != None:
        match = re.search(bus_regex, params).group(0)
        read_id = int(match.split(":")[-1].split(",")[0])
    elif re.search(bus_regex_short, params) != None:
        tech = re.search(bus_regex_short, params).group(0)
        read_id = kb_tech_dict[tech.lower()]
    else:
        logger.error("Not a valid BUS(barcode:umi:set) string. Please check -x argument")
        os._exit(1)  # noqa
    return read_id


def hex_to_rgb(value):
    """In: hex(#ffffff). Out: tuple(255, 255, 255)"""
    value = value.lstrip("#")
    rgb = tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))
    return rgb


def rgb_to_hsv(value):
    """In: tuple(1, 1, 0.996) or tuple(255, 255, 254). Out: tuple(0.24, 0.04, 1)"""
    if not all([n <= 1 for n in value]):
        value = tuple([n / 255 for n in value])
    hsv = colorsys.rgb_to_hsv(value[0], value[1], value[2])
    return hsv


def hsv_to_ucsc(value):
    """
    UCSC accepts RGB as string without spaces.
    In: tuple(1, 1, 0.996). Out: str("255,255,254")
    """
    # older versions of numpy hijack round and return a float, hence int()
    # see https://github.com/numpy/numpy/issues/11810
    rgb = [int(round(n * 255)) for n in mcolors.hsv_to_rgb(value)]
    ucsc_rgb = f"{rgb[0]},{rgb[1]},{rgb[2]}"
    return ucsc_rgb


def color_parser(color: str, color_dicts: list = None) -> tuple:
    """
    convert a string with RGB/matplotlib named colors to matplotlib HSV tuples.

    supports RGB colors with ranges between 0-1 or 0-255.

    supported matplotlib colors can be found here:
    https://matplotlib.org/3.3.1/gallery/color/named_colors.html
    """
    # input: RGB
    if color.count(",") == 2:
        value = [float(c) for c in color.split(",")]
        return rgb_to_hsv(value)

    # input: matplotlib colors
    cdicts = color_dicts if color_dicts else DEFAULT_COLOR_DICTS
    for cdict in cdicts:
        if color in cdict:
            value = cdict[color]

            # tableau, css4 and xkcd return hex colors.
            if str(value).startswith("#"):
                value = hex_to_rgb(value)

            return rgb_to_hsv(value)

    logger.error(f"Color not recognized: {color}")
    os._exit(1)  # noqa


def color_picker(n, min_h=0, max_h=0.85, s=1.00, v=0.75, alternate=True):
    """
    Return a list of n tuples with HSV colors varying only in hue.
    "Alternate" determines whether hues transition gradually or discretely.
    """
    # for fewer samples, select nearby colors
    steps = max(n, 8)

    hues = np.linspace(min_h, max_h, steps).tolist()[0:n]
    if alternate:
        m = ceil(len(hues) / 2)
        h1 = hues[:m]
        h2 = hues[m:]
        hues[::2] = h1
        hues[1::2] = h2

    hsv_colors_list = [(h, s, v) for h in hues]
    return hsv_colors_list


def color_gradient(hsv: tuple, n: int) -> List[tuple]:
    """
    Based on the input HSV color,
    return a list of n tuples with HSV colors
    of increasing brightness (saturation + value).
    """
    # for fewer samples, select nearby colors
    steps = max(n, 4)

    h = hsv[0]
    s = np.linspace(hsv[1], 0.2, steps)  # goes down
    v = np.linspace(hsv[2], 1.0, steps)  # goes up

    hsv_gradient_list = [(h, s[i], v[i]) for i in range(n)]
    return hsv_gradient_list


def shorten(string, max_length, methods="right"):
    """
    shorten a string to a max_length, multiple methods accepted.
    "signs" and "vowels" remove their respective characters from right to left.
    "left","right" and "center" can be performed afterward if the desired length is not yet reached.
    """
    overhead = len(string) - max_length
    if overhead <= 0:
        return string

    if "signs" in methods:
        signs = ["-", "_", "."]
        s = ""
        for l in string[::-1]:
            if l in signs and overhead > 0:
                overhead -= 1
            else:
                s += l
        string = s[::-1]

    if "vowels" in methods:
        vowels = ["a", "e", "i", "o", "u"]
        s = ""
        for l in string[::-1]:
            if l in vowels and overhead > 0:
                overhead -= 1
            else:
                s += l
        string = s[::-1]

    if "right" in methods:
        string = string[:max_length]
    elif "left" in methods:
        string = string[len(string) - max_length :]
    elif "center" in methods:
        string = string[: ceil(max_length / 2)] + string[len(string) - floor(max_length / 2) :]

    return string


class CaptureStdout(list):
    """
    Context manager that somehow manages to capture prints,
    and not snakemake log
    """

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


class CaptureStderr(list):
    """
    Context manager that somehow manages to capture prints,
    and not snakemake log
    """

    def __enter__(self):
        self._stderr = sys.stderr
        sys.stderr = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stderr = self._stderr


def prep_filelock(lock_file, max_age=10):
    """
    create the directory for the lock_file if needed
    and remove locks older than the max_age (in seconds)
    """
    os.makedirs(os.path.dirname(lock_file), exist_ok=True)

    # sometimes two jobs start in parallel and try to delete at the same time
    try:
        # ignore locks that are older than the max_age
        if os.path.exists(lock_file) and time.time() - os.stat(lock_file).st_mtime > max_age:
            os.unlink(lock_file)
    except FileNotFoundError:
        pass


def retry_pickling(func):
    def wrap(*args, **kwargs):
        # we get two tries, in case parallel executions are interfering with one another
        for _ in range(2):
            try:
                return func(*args, **kwargs)
            except FileNotFoundError:
                time.sleep(1)
        else:
            logger.error("There were some problems with locking the seq2science cache. Please try again in a bit.")
            os._exit(1)  # noqa

    return wrap


class PickleDict(dict):
    """dict with builtin pickling utility"""

    def __init__(self, file):
        self.file = file
        self.filelock = f"{file}.p.lock"

        data = self.load() if os.path.exists(self.file) else dict()
        super(PickleDict, self).__init__(data)

    @retry_pickling
    def load(self):
        prep_filelock(self.filelock, 30)
        with FileLock(self.filelock):
            return pickle.load(open(self.file, "rb"))

    @retry_pickling
    def save(self):
        prep_filelock(self.filelock, 30)
        with FileLock(self.filelock):
            pickle.dump(self, open(self.file, "wb"))

    def search(self, search_assemblies: list):
        """
        Check the genomepy database for a provider stat serves both genome and annotation for an assembly.
        If impossible, settle with the first provider that serves the genome.
        """
        logger.info("Determining assembly providers, this can take some time.")
        for assembly in search_assemblies:
            if assembly not in self:
                self[assembly] = {"genome": None, "annotation": None}

        with open(logger.logfile, "w") as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                genomepy.logger.remove()
                genomepy.logger.add(
                    log,
                    format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
                    level="INFO",
                )
                for p in genomepy.providers.online_providers():
                    search_assemblies = [a for a in search_assemblies if self[a]["annotation"] is None]
                    for assembly in search_assemblies:
                        if assembly not in p.genomes:
                            continue  # check again next provider

                        if p.annotation_links(assembly):
                            self[assembly]["genome"] = p.name
                            self[assembly]["annotation"] = p.name

                        elif self[assembly]["genome"] is None:
                            self[assembly]["genome"] = p.name

                    if all(self[a]["annotation"] for a in search_assemblies):
                        break  # don't load the next provider

        # store added assemblies
        self.save()

    def check(self, assemblies: list, annotation_required: bool, verbose: bool):
        """
        Check if the genome (and the annotation if required) can be downloaded for each specified assembly.
        """
        for assembly in assemblies:
            if self[assembly]["genome"] is None:
                logger.warning(
                    f"Could not download assembly {assembly}.\n"
                    f"Find alternative assemblies with `genomepy search {assembly}`"
                )
                os._exit(1)  # noqa

            if self[assembly]["annotation"] is None:
                if verbose:
                    logger.warning(
                        f"No annotation for assembly {assembly} can be downloaded. Another provider (and "
                        f"thus another assembly name) might have a gene annotation.\n"
                        f"Find alternative assemblies with `genomepy search {assembly}`\n"
                    )
                if annotation_required:
                    os._exit(1)  # noqa


def is_local(assembly: str, ftype: str, config: dict) -> bool:
    """checks if genomic file(s) are present locally"""
    file = os.path.join(config["genome_dir"], assembly, assembly)
    local_fasta = os.path.exists(f"{file}.fa")
    local_gtf = os.path.exists(f"{file}.annotation.gtf")
    local_bed = os.path.exists(f"{file}.annotation.bed")
    if ftype == "genome":
        return local_fasta
    if ftype == "annotation":
        # check genome and annotations, as genome is always needed
        return local_gtf and local_bed and local_fasta


def _get_yaml_file(rules_dir):
    """
    Return the path to the requirements.yaml file
    If s2s is installed in editable mode, the editable yaml file is returned.
    """
    # first check if installed normally
    yaml_file = os.path.abspath(os.path.join(
        rules_dir, "..", "envs", "seq2science_requirements.yaml"
    ))
    if os.path.isfile(yaml_file):
        return yaml_file

    # otherwise check if installed in editable mode
    yaml_file = os.path.abspath(os.path.join(
        rules_dir, "..", "..", "requirements.yaml"
    ))
    if os.path.isfile(yaml_file):
        return yaml_file

    raise FileNotFoundError(
        "Seq2science couldn't find it's own requirements file! "
        "This shouldn't happen, so please raise an issue on github!"
    )


def _get_yaml_versions(yaml_file):
    """
    Return a dict with packages and versions from requirements.yaml.
    """
    with open(yaml_file, "r") as stream:
        env = yaml.safe_load(stream)

    versions = {}
    for dependency in env["dependencies"]:
        # remove channel prefix
        if "::" in dependency:
            dependency = dependency.split("::")[1]
        # split tool and version (ignore build if present)
        package, version = dependency.split("=")[0:2]
        versions[package] = version
    return versions


def _get_current_version(package):
    """
    Attempt to return a given package's version
    """
    # package-isolation is not a package
    # xdg keeps its version in a pyproject.toml (not included)
    # argcomplete keeps its version in a setup.py (not included)
    # trackhub versioning is weird
    if package in ["conda-ecosystem-user-package-isolation", "xdg", "argcomplete", "trackhub"]:
        return None
    if package == "python":
        return sys.version.split()[0]

    # some packages have different names on conda
    if package == "snakemake-minimal":
        package = "snakemake"
    elif package == "pyyaml":
        package = "yaml"
    elif package == "biopython":
        package = "Bio"

    ldict = dict()
    exec(f"from {package} import __version__", {}, ldict)
    current_version = ldict["__version__"]
    return current_version


def assert_versions(rules_dir):
    """
    For each package, check that the installed version matches the required version
    """
    error = False
    yaml_file = _get_yaml_file(rules_dir)
    versions = _get_yaml_versions(yaml_file)
    for package, required_version in versions.items():
        current_version = _get_current_version(package)
        if current_version is None:
            continue
        if not current_version.startswith(required_version):
            logger.error(
                f"Seq2science requires {package.capitalize()} version {required_version}, "
                f"found version {current_version}."
            )
            error = True
    if error:
        logger.error("Please create a new conda environment.\n")
        os._exit(1)  # noqa
