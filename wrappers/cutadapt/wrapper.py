"""Snakemake wrapper for trimming paired-end reads using cutadapt."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
from detect_adapter import detect_adapters_and_cnts, adapters
import gzip
import re


def is_interleaved(fname):
    NLINES = 40
    with gzip.open(fname) as f:
        ids = [f.readline() for i in range(NLINES)]
    ids = [re.split('/| ', ids[i].decode())[0] for i in range(0, NLINES, 4)]

    return len(ids) / 2 == len(set(ids))

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if hasattr(snakemake.params, "adapter"):
    if snakemake.params.adapter == "auto":
        result = detect_adapters_and_cnts(snakemake.input[0])
        if len(result[0]) > 0:
            snakemake.params.adapter = adapters[result[0][0]].decode()
        else: 
            snakemake.params.adapter = ""
    if  snakemake.params.adapter != "":
        snakemake.params.extra += " -a {}".format(snakemake.params.adapter)

if is_interleaved(snakemake.input[0]):
    snakemake.params.extra += " --interleaved"

shell(
    "cutadapt"
    " {snakemake.params.extra}"
    " -j {snakemake.threads}"   
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}")
