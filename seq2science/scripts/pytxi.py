"""
convert TPM to gene counts
"""
import os
import contextlib
import pytxi


logfile = snakemake.log[0]
fnames = snakemake.input.cts
assembly = snakemake.wildcards.assembly
species = snakemake.input.fa[0]
output = snakemake.output[0]
sample_names = snakemake.params.names

# redirect all messages to a logfile
with open(logfile, "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        pytxi.logger.remove()
        pytxi.logger.add(logfile)

        fnames = [f"{d}/quant.sf" for d in fnames]
        outdir = os.path.dirname(output)
        os.makedirs(outdir, exist_ok=True)

        txi = pytxi.TxImport()
        txi.import_files(fnames, sample_names=sample_names, tx2gene=None, species=species)
        txi.abundance.to_csv(f"{outdir}/{assembly}-TPM.tsv", sep="\t")
        txi.counts.to_csv(f"{outdir}/{assembly}-counts.tsv", sep="\t")
        txi.length.to_csv(f"{outdir}/{assembly}-gene_length.tsv", sep="\t")
