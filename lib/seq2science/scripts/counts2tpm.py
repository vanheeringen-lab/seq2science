import contextlib

import numpy as np
import pandas as pd
from genomepy import Annotation


counts_file = snakemake.input.cts[0]
gtf = snakemake.input.gtf[0]
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
        ann = Annotation(gtf)
        counts = pd.read_table(counts_file, comment="#", index_col=0)

        # what are the counts called: gene names or gene ids?
        counts_names = set(counts.index)
        test = ann.gtf["attribute"].head(100)
        name_overlap = 0
        if any(test.str.contains("gene_name")):
            gtf_names = set(ann.from_attributes("gene_name"))
            name_overlap = len(counts_names & gtf_names)
        id_overlap = 0
        if any(test.str.contains("gene_id")):
            gtf_names = set(ann.from_attributes("gene_id"))
            id_overlap = len(counts_names & gtf_names)
        if name_overlap == 0 and id_overlap == 0:
            raise ValueError(
                "names in the counts table don't correspond to either the GTF's "
                "gene_name or gene_id!"
            )
        print(
            f"{int(100*max(name_overlap,id_overlap)/len(counts_names))}% of names "
            "overlap between counts and GTF."
        )
        attribute = "gene_name" if name_overlap > id_overlap else "gene_id"

        lengths = ann.lengths(attribute)
        tpms = counts2tpm(counts.join(lengths, how="inner"))

        lengths.to_csv(out_lengths, sep="\t")
        tpms.to_csv(out_tpms, sep="\t")
