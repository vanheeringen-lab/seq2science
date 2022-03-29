"""
convert TPM to gene counts
"""
import os
import contextlib

from genomepy import Annotation, Genome
import pandas as pd
import pytxi


logfile = snakemake.log[0]
fnames = snakemake.input.cts
assembly = snakemake.wildcards.assembly
fa = snakemake.input.fa[0]
gtf = snakemake.input.gtf[0]
from_gtf = snakemake.params.from_gtf
out_counts = snakemake.output.counts[0]
out_tpms = snakemake.output.tpms[0]
out_lengths = snakemake.output.lengths[0]
sample_names = snakemake.params.names

# redirect all messages to a logfile
open(logfile, "w")  # start a clean log
with open(logfile, "a") as log:  # appending because we mix streams (tqdm, logger & std)
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        pytxi.logger.remove()
        pytxi.logger.add(logfile)

        fnames = [f"{d}/quant.sf" for d in fnames]
        outdir = os.path.dirname(out_counts)
        os.makedirs(outdir, exist_ok=True)

        tx2gene = None
        if from_gtf:
            ann = Annotation(gtf)

            # check if we can convert transcript ids to symbols with the GTF file
            test = ann.gtf["attribute"].head(100)
            attribute = None
            if any(test.str.contains("gene_name")):
                attribute = "gene_name"
            elif any(test.str.contains("gene_id")):
                attribute = "gene_id"

            if attribute is not None:
                dct = ann.gtf_dict("transcript_id", attribute)
                df = pd.DataFrame.from_dict(dct, orient="index", columns=[attribute])
                df.index.name = "transcript_id"
                tx2gene = f"{outdir}/tx2gene.tsv"
                df.to_csv(tx2gene, sep="\t")

                # print some metrics
                lendf = len(set(df.index))
                lenbed = len(ann.genes("bed"))
                pytxi.logger.info(
                    f"GTF has {lendf} transcript ids, TPMs have {lenbed} transcript ids"
                )
            else:
                print("GTF does not contain gene_name or gene_id! Using MyGene.info...")

        taxonomy = None
        if tx2gene is None:
            g = Genome(fa, build_index=False)
            taxonomy = g.tax_id
            if not isinstance(taxonomy, int):
                raise ValueError(f"{g.name} README.txt does not contain a valid taxonomy id!")

        txi = pytxi.TxImport()
        txi.import_files(fnames, sample_names, tx2gene, taxonomy)

        txi.abundance.to_csv(out_tpms, sep="\t")
        txi.counts.to_csv(out_counts, sep="\t")
        txi.length.to_csv(out_lengths, sep="\t")
