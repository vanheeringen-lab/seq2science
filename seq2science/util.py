"""
Utility functions for seq2science
"""
import colorsys
import os
import sys
from typing import List
from math import ceil, floor
import time
import urllib.request
import numpy as np
import pandas as pd
import pysradb
import matplotlib.colors as mcolors
from snakemake.io import expand


def _sample_to_idxs(df: pd.DataFrame, sample: str) -> List[int]:
    """
    Get a list of index/indices that belong to a run in a pysradb dataframe
    """
    if sample.startswith(("SRR", "DRR", "ERR")):
        idxs = df.index[df.run_accession == sample].tolist()
        assert len(idxs) == 1, f"sample {sample} with idxs: {idxs}"
    elif sample.startswith("SRX"):
        idxs = df.index[df.experiment_accession == sample].tolist()
        assert len(idxs) >= 1, len(idxs)
    else:
        assert False, f"sample {sample} not a run, this should not be able to happen!" \
                      f" Please make an issue about this!"
    return idxs


def samples2metadata_local(samples: List[str], config: dict) -> dict:
    """

    """
    sampledict = dict()
    for sample in samples:
        if os.path.exists(expand(f'{{fastq_dir}}/{sample}.{{fqsuffix}}.gz', **config)[0]):
            sampledict[sample] = dict()
            sampledict[sample]["layout"] = "SINGLE"
        elif all(os.path.exists(path) for path in expand(f'{{fastq_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz', **config)):
            sampledict[sample] = dict()
            sampledict[sample]["layout"] = "PAIRED"
        elif sample.startswith(('GSM', 'SRX', 'SRR', 'ERR', 'DRR')):
            continue
        else:
            raise ValueError(f"\nsample {sample} was not found..\n"
                             f"We checked for SE file:\n"
                             f"\t{config['fastq_dir']}/{sample}.{config['fqsuffix']}.gz \n"
                             f"and for PE files:\n"
                             f"\t{config['fastq_dir']}/{sample}_{config['fqext1']}.{config['fqsuffix']}.gz \n"
                             f"\t{config['fastq_dir']}/{sample}_{config['fqext2']}.{config['fqsuffix']}.gz \n"
                             f"and since the sample did not start with either GSM, SRX, SRR, ERR, and DRR we "
                             f"couldn't find it online..\n")
    return sampledict


def samples2metadata_sra(samples: List[str]) -> dict:
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
    sampledict = {sample: dict() for sample in samples}

    # only continue with public samples
    db_sra = pysradb.SRAweb()

    # all samples that are on GEO (GSM numbers), must first be converted to SRA ids (SRX numbers)
    geo_samples = [sample for sample in samples if sample.startswith("GSM")]

    # in sample2clean we store the (potential GEO) sample name in a SRA compliant name
    try:
        df_geo = db_sra.gsm_to_srx(geo_samples)
    except:
        print("We had trouble querying the SRA. This probably means that the SRA was unresponsive, and their servers "
              "are overloaded or slow. Please try again in a bit...")
        sys.exit(1)

    sample2clean = dict(zip(df_geo.experiment_alias, df_geo.experiment_accession))

    # now add the already SRA compliant names with a reference to itself
    sample2clean.update({sample: sample for sample in samples if sample not in geo_samples})

    # check our samples on sra
    try:
        df_sra = db_sra.sra_metadata(list(sample2clean.values()), detailed=True)
    except:
        print("We had trouble querying the SRA. This probably means that the SRA was unresponsive, and their servers "
              "are overloaded or slow. Please try again in a bit...")
        sys.exit(1)

    for sample, clean in sample2clean.items():
        # table indices
        idxs = _sample_to_idxs(df_sra, clean)

        # get all runs that belong to the sample
        runs = df_sra.loc[idxs].run_accession.tolist()
        assert len(runs) >= 1
        sampledict[sample]["runs"] = runs

        # get the layout
        layout = df_sra.loc[idxs].library_layout.tolist()
        assert len(set(layout)) == 1, f"sample {sample} consists of mixed layouts, bad!"
        assert layout[0] in ["PAIRED", "SINGLE"], f"sample {sample} is an unclear layout, bad!"
        sampledict[sample]["layout"] = layout[0]

        # get the ena url
        sampledict[sample]["ena_fastq_ftp"] = dict()
        for run in runs:
            if layout[0] == "SINGLE":
                sampledict[sample]["ena_fastq_ftp"][run] = df_sra[df_sra.run_accession == run].ena_fastq_ftp.tolist()
            elif layout[0] == "PAIRED":
                sampledict[sample]["ena_fastq_ftp"][run] = df_sra[df_sra.run_accession == run].ena_fastq_ftp_1.tolist() + df_sra[
                    df_sra.run_accession == run].ena_fastq_ftp_2.tolist()

        # if any run from a sample is not found on ENA, better be safe, and assume that sample as a whole is not on ENA
        if any(["N/A" in urls for run, urls in sampledict[sample]["ena_fastq_ftp"].items()]):
            sampledict[sample]["ena_fastq_ftp"] = None

    return sampledict


def samples2metadata(samples: List[str], config: dict) -> dict:
    local_samples = samples2metadata_local(samples, config)
    public_samples = [sample for sample in samples if sample not in local_samples.keys()]

    if len(public_samples) == 0:
        return local_samples

    sra_samples = samples2metadata_sra(public_samples)

    return {**local_samples, **sra_samples}


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


def sieve_bam(configdict):
    """
    helper function to check whether or not we use rule sieve_bam
    """
    return (
        configdict.get("min_mapping_quality", 0) > 0
        or configdict.get("tn5_shift", False)
        or configdict.get("remove_blacklist", False)
        or configdict.get("remove_mito", False)
    )


def parse_de_contrasts(de_contrast):
    """
    Extract contrast column, groups and batch effect column from a DE contrast design

    Accepts a string containing a DESeq2 contrast design

    Returns
    parsed_contrast, lst: the contrast components
    batch, str: the batch effect column or None
    """
    # remove whitespaces (and '~'s if used)
    de_contrast = de_contrast.replace(" ", "").replace("~", "")

    # split and store batch effect
    batch = None
    if "+" in de_contrast:
        batch, de_contrast = de_contrast.split("+")[0:2]

    # parse contrast column and groups
    parsed_contrast = de_contrast.split("_")
    return parsed_contrast, batch


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


def hex_to_rgb(value):
    """hex in, tuple out (0-255)"""
    value = value.lstrip('#')
    return tuple(int(value[i:i + 2], 16) for i in (0, 2, 4))


def rgb_to_hvs(value):
    """In: tuple(1, 1, 0.996) or tuple(255, 255, 254). Out: tuple(0.24, 0.04, 1)"""
    if not all([n <= 1 for n in value]):
        value = tuple([n / 255 for n in value])
    return colorsys.rgb_to_hsv(value[0], value[1], value[2])


def hsv_to_ucsc(value):
    """
    UCSC accepts RGB as string without spaces.
    In: tuple(1, 1, 0.996). Out: str("255,255,254")
    """
    rgb = [round(n*255) for n in mcolors.hsv_to_rgb(value)]
    return f"{rgb[0]},{rgb[1]},{rgb[2]}"


# default colors in matplotlib. Order dictates the priority.
default_color_dicts = [mcolors.BASE_COLORS, mcolors.TABLEAU_COLORS, mcolors.CSS4_COLORS, mcolors.XKCD_COLORS]


def color_parser(color, color_dicts=None):
    """
    convert RGB/matplotlib named colors to matplotlib HSV tuples.

    supported matplotlib colors can be found here:
    https://matplotlib.org/3.3.1/gallery/color/named_colors.html
    """
    # input: RGB
    if color.count(",") == 2:
        value = [float(c) for c in color.split(",")]
        return rgb_to_hvs(value)

    # input: matplotlib colors
    cdicts = color_dicts if color_dicts else default_color_dicts
    for cdict in cdicts:
        if color in cdict:
            value = cdict[color]

            # tableau, css4 and xkcd return hex colors.
            if str(value).startswith("#"):
                value = hex_to_rgb(value)

            return rgb_to_hvs(value)

    print(f"Color not recognized: {color}")
    raise ValueError


def color_picker(n, min_h=0, max_h=0.85, dflt_s=1.00, dflt_v=0.75, alternate=True):
    """
    Return a list of n tuples with HSV colors varying only in hue.
    "Alternate" determines whether hues transition gradually or discretely.
    """
    # for fewer samples, select nearby colors
    steps = max(n, 8)

    hues = np.linspace(min_h,max_h, steps).tolist()[0:n]
    if alternate:
        m = ceil(len(hues)/2)
        h1 = hues[:m]
        h2 = hues[m:]
        hues[::2] = h1
        hues[1::2] = h2

    S = dflt_s
    V = dflt_v
    return [(H, S, V) for H in hues]


def color_gradient(hsv, n):
    """
    Based on the input HSV color,
    return a list of n tuples with HSV colors varying in brightness (saturation + value).
    """
    # for fewer samples, select nearby colors
    steps = max(n, 4)

    H = [hsv[0] for _ in range(n)]       # constant
    S = np.linspace(hsv[1], 0.2, steps)  # goes down
    V = np.linspace(hsv[2], 1.0, steps)  # goes up
    return [(H[_], S[_], V[_]) for _ in range(n)]


def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]


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
        string = string[len(string)-max_length:]
    elif "center" in methods:
        string = string[:ceil(max_length/2)] + string[len(string)-floor(max_length/2):]

    return string
