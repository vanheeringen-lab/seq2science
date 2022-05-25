from contextlib import redirect_stdout, redirect_stderr
from os.path import dirname, join
from shutil import copyfile

# human and mouse are supported by the default m2f
supported_genome_prefixes = ["GRCh", "GRCm", "hg19", "hg38", "mm10", "mm39"]

# log errors
with open(str(snakemake.log), "w") as f:
    with redirect_stdout(f), redirect_stderr(f):

        if snakemake.input.genome[:4] in supported_genome_prefixes:
            # copy the default m2f
            from gimmemotifs.motif import pfmfile_location

            in_pfmfile = pfmfile_location(None)
            out_pfmfile = snakemake.output[0]

            in_m2ffile = in_pfmfile.replace(".pfm", ".motif2factors.txt")
            out_m2ffile = out_pfmfile.replace(".pfm", ".motif2factors.txt")

            copyfile(in_pfmfile, out_pfmfile)
            copyfile(in_m2ffile, out_m2ffile)

        else:
            # create an ortholog m2f
            from gimmemotifs.orthologs import motif2factor_from_orthologs

            motif2factor_from_orthologs(
                database=snakemake.params.database,
                new_reference=[snakemake.input.genome],
                extra_orthologs_references=[snakemake.params.motif2factors_reference],
                genomes_dir=snakemake.params.genomes_dir,
                outdir=dirname(snakemake.output[0]),
                threads=snakemake.threads,
            )
