"""
convert TPM to gene counts
"""
import os
import contextlib

from genomepy import Annotation
import pandas as pd
import pytxi


logfile = snakemake.log[0]
fnames = snakemake.input.cts
assembly = snakemake.wildcards.assembly
species = snakemake.input.fa[0]
from_gtf = snakemake.params.from_gtf
output = snakemake.output[0]
sample_names = snakemake.params.names

# redirect all messages to a logfile
open(logfile, "w")  # start a clean log
with open(logfile, "a") as log:  # appending because we mix streams (tqdm, logger & std)
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        pytxi.logger.remove()
        pytxi.logger.add(logfile)

        fnames = [f"{d}/quant.sf" for d in fnames]
        outdir = os.path.dirname(output)
        os.makedirs(outdir, exist_ok=True)

        tx2gene = None
        if from_gtf:
            # check if we can convert transcript ids to symbols (gene names) with the GTF file
            ann = Annotation(species)
            dct = ann.gtf_dict("transcript_id", "gene_name")
            df = pd.DataFrame.from_dict(dct, orient="index", columns=["gene_name"])
            df.index.name = "transcript_id"
            lendf = len(set(df.index))
            lengenes = len(ann.genes("bed"))
            pytxi.logger.info(
                f"GTF has {lendf} transcript ids, TPMs have {lengenes} transcript ids"
            )
            if lendf >= 0.9*lengenes:
                tx2gene = f"{outdir}/tx2gene.tsv"
                df.to_csv(tx2gene, sep="\t")

        txi = pytxi.TxImport()
        txi.import_files(fnames, sample_names, tx2gene, species)
        txi.abundance.to_csv(f"{outdir}/{assembly}-TPM.tsv", sep="\t")
        txi.counts.to_csv(f"{outdir}/{assembly}-counts.tsv", sep="\t")
        txi.length.to_csv(f"{outdir}/{assembly}-gene_length.tsv", sep="\t")
