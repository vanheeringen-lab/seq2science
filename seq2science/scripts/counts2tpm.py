import contextlib

import numpy as np
import pandas as pd
from genomepy import Annotation


counts_file = snakemake.input.cts[0]
species = snakemake.input.gtf[0]
out_tpms = snakemake.output.tpms[0]
out_lengths = snakemake.output.lengths[0]
log_file = snakemake.log[0]


def counts2tpm(df):
    """
    convert read counts to TPM (transcripts per million)
    """
    gene_names = df.index
    sample_names = [c for c in df.columns if c != "length"]

    sample_reads = df[sample_names]
    gene_length = df[["length"]]
    rate = sample_reads.values / gene_length.values
    tpm = rate / np.sum(rate, axis=0).reshape(1, -1) * 1e6
    return pd.DataFrame(data=tpm, columns=sample_names, index=gene_names)


# redirect all messages to a logfile
open(log_file, "w")  # start a clean log
with open(log_file, "a") as log:  # appending because we mix streams (tqdm, logger & std)
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        # check if we can convert transcript ids to symbols (gene names) with the GTF file
        ann = Annotation(species)

        counts = pd.read_table(counts_file, comment="#", index_col=0)
        lengths = ann.lengths("gene_name")
        tpms = counts2tpm(counts.join(lengths, how="inner"))

        lengths.to_csv(out_tpms, sep="\t")
        tpms.to_csv(out_tpms, sep="\t")
