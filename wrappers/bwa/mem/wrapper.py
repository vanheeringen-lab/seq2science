__author__ = "Johannes Köster, Julian de Ruiter"
__copyright__ = "Copyright 2016, Johannes Köster and Julian de Ruiter"
__email__ = "koester@jimmy.harvard.edu, julianderuiter@gmail.com"
__license__ = "MIT"


from os import path
import gzip
from snakemake.shell import shell
import re


def is_interleaved(fname):
    NLINES = 40
    with gzip.open(fname) as f:
        ids = [f.readline() for i in range(NLINES)]
    ids = [re.split('/| ', ids[i].decode())[0] for i in range(0, NLINES, 4)]

    return len(ids) / 2 == len(set(ids))


# Extract arguments.
extra = snakemake.params.get("extra", "")

sort = snakemake.params.get("sort", "none")
sort_order = snakemake.params.get("sort_order", "coordinate")
sort_extra = snakemake.params.get("sort_extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Check inputs/arguments.
if len(snakemake.input.reads) not in {1, 2}:
    raise ValueError("input must have 1 (single-end) or "
                     "2 (paired-end) elements")

if sort_order not in {"coordinate", "queryname"}:
    raise ValueError("Unexpected value for sort_order ({})".format(sort_order))

# Determine which pipe command to use for converting to bam or sorting.
if sort == "none":

    # Simply convert to bam using samtools view.
    pipe_cmd = "samtools view -Sbh -o {snakemake.output[0]}"

elif sort == "sambamba":

    # Sort alignments using sambamba sort.
    threads = sorted([1, 5, snakemake.threads // 3])[1]
    pipe_cmd = f"sambamba view --nthreads {threads} -S -f bam -o /dev/stdout /dev/stdin | " \
        f"sambamba sort /dev/stdin {{sort_extra}} -o {{snakemake.output[0]}} --nthreads {threads}"

    # Add name flag if needed.
    if sort_order == "queryname":
        sort_extra += " -n"

elif sort == "samtools":

    # Sort alignments using samtools sort.
    threads = sorted([1, 4, snakemake.threads // 3])[1]
    pipe_cmd = f"samtools sort {{sort_extra}} -o {{snakemake.output[0]}} -m 500M --threads {threads}"

    # Add name flag if needed.
    if sort_order == "queryname":
        sort_extra += " -n"

    prefix = path.splitext(snakemake.output[0])[0]
    sort_extra += " -T " + prefix + ".tmp"

elif sort == "picard":
    assert False, "picard not implemented in environment"
    # Sort alignments using picard SortSam.
    pipe_cmd = ("picard SortSam {sort_extra} INPUT=/dev/stdin"
                " OUTPUT={snakemake.output[0]} SORT_ORDER={sort_order}")

else:
    raise ValueError("Unexpected value for params.sort ({})".format(sort))

if is_interleaved(snakemake.input.reads[0]):
    extra += " -p"

# format the remaining wildcards
snakemake.params.index = snakemake.params.index[0].format(**dict(snakemake.wildcards))

shell(
    f"exec {{log}}; "
    f"trap \"rm -f {'.'.join(snakemake.output[0].split('.')[:-1])}*\" INT;"
    f"cpulimit --include-children -l {{snakemake.threads}}00 -- "
    f"bwa mem"
    f" -t {{snakemake.threads}}"
    f" {{extra}}"
    f" {{snakemake.params.index}}"
    f" {{snakemake.input.reads}}"
    f" | " + pipe_cmd)
