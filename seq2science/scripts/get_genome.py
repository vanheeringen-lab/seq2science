"""
Script to download genome
"""
import contextlib

import genomepy
from snakemake.logging import logger


with open(snakemake.log[0], "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        # select a provider with the annotation if possible
        a = snakemake.params.providers[snakemake.wildcards.raw_assembly]["annotation"]
        g = snakemake.params.providers[snakemake.wildcards.raw_assembly]["genome"]
        provider = g if a is None else a
        if provider not in ["ensembl", "ucsc", "ncbi"]:
            logger.error("Seq2science database error. Please run 'seq2science --clean' and try again")
            raise TerminatedException

        p = genomepy.ProviderBase.create(provider)
        p.download_genome(snakemake.wildcards.raw_assembly, snakemake.params.genome_dir)

        # try to download the blacklist
        genome = genomepy.Genome(snakemake.wildcards.raw_assembly, snakemake.params.genome_dir)
        plugins = genomepy.plugin.init_plugins()
        plugins["blacklist"].after_genome_download(genome)
