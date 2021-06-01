"""
Script to download genome
"""
import os
import contextlib

import genomepy


with open(snakemake.log[0], "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        # list user plugins
        active_plugins = genomepy.functions.config.get("plugin", [])
        # deactivate user plugins
        genomepy.functions.manage_plugins("disable", active_plugins)

        # select a provider with the annotation if possible
        a = snakemake.params.providers[snakemake.wildcards.raw_assembly]["annotation"]
        g = snakemake.params.providers[snakemake.wildcards.raw_assembly]["genome"]
        provider = g if a is None else a
        try:
            genomepy.install_genome(
                name=snakemake.wildcards.raw_assembly,
                provider=provider,
                genomes_dir=snakemake.params.genome_dir,
                force=True,
            )

            # delete the support files
            # (we recreate these separately to make the output simple for snakemake and prevent redownloading)
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

        finally:
            # reactivate user plugins
            genomepy.functions.manage_plugins("enable", active_plugins)
