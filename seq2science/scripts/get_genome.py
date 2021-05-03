"""
Script to download genome
"""
import os
import contextlib

import genomepy


with open(snakemake.log[0], "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        try:
            # select a provider with the annotation if possible
            a = snakemake.params.providers[snakemake.wildcards.raw_assembly]["annotation"]
            g = snakemake.params.providers[snakemake.wildcards.raw_assembly]["genome"]
            provider = g if a is None else a

            p = genomepy.ProviderBase.create(provider)
            p.download_genome(snakemake.wildcards.raw_assembly, snakemake.params.genome_dir)

            # try to download the blacklist
            genome = genomepy.Genome(snakemake.wildcards.raw_assembly, snakemake.params.genome_dir)
            plugins = genomepy.plugin.init_plugins()
            plugins["blacklist"].after_genome_download(genome)

            # now delete the support files as these are created by a different rule
            os.remove(f"{snakemake.output[0]}.fai")
            os.remove(f"{snakemake.output[0]}.sizes")
            os.remove(f"{snakemake.output[0][:-2]}gaps.bed")
        except Exception as e:
            print(e)
            print("\nSomething went wrong while downloading the genome (see error message above). "
                  "When this happens it is almost always because we had troubles connecting to the"
                  "servers hosting the genome assemblies. Usually this is resolved by just running seq2science"
                  "again, either immediately or in a couple hours.\n\n"
                  "If the problem persists you could try running `seq2science clean` and see if that resolves the "
                  "issue.")
