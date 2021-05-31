"""
Script to download the genome blacklist.

If none exists, an empty file is generated to keep snakemake content.
"""
import os
import contextlib

import genomepy


with open(snakemake.log[0], "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        genome = genomepy.Genome(snakemake.wildcards.raw_assembly, snakemake.params.genome_dir)
        plugins = genomepy.plugin.init_plugins()
        plugins["blacklist"].after_genome_download(genome)

        # touch the file if no blacklist was available
        if not os.path.exists(snakemake.output[0]):
            with open(snakemake.output[0], 'a'):
                os.utime(snakemake.output[0], None)
        #
        # # list user plugins
        # active_plugins = genomepy.functions.config.get("plugin", [])
        # # deactivate user plugins
        # genomepy.funcions.manage_plugins("disable", active_plugins)
        # # activate rule plugin (blacklist)
        # genomepy.funcions.manage_plugins("enable", ["blacklist"])
        #
        # # select a provider with the annotation if possible
        # a = snakemake.params.providers[snakemake.wildcards.raw_assembly]["annotation"]
        # g = snakemake.params.providers[snakemake.wildcards.raw_assembly]["genome"]
        # provider = g if a is None else a
        # try:
        #     genomepy.install_genome(
        #         name=snakemake.wildcards.raw_assembly,
        #         provider=provider,
        #         genomes_dir=snakemake.params.genome_dir,
        #         force=False,
        #     )
        #
        #     # touch the file if no blacklist was available
        #     if not os.path.exists(snakemake.output[0]):
        #         with open(snakemake.output[0], 'a'):
        #             os.utime(snakemake.output[0], None)
        # except Exception as e:
        #     print(e)
        #     print("\nSomething went wrong while downloading the blacklist (see error message above). "
        #           "When this happens it is almost always because we had troubles connecting to the"
        #           "servers hosting the genome assemblies. Usually this is resolved by just running seq2science"
        #           "again, either immediately or in a couple hours.\n\n"
        #           "If the problem persists you could try running `seq2science clean` and see if that resolves the "
        #           "issue.")
        #
        # finally:
        #     # deactivate rule plugin (blacklist)
        #     genomepy.funcions.manage_plugins("disable", ["blacklist"])
        #     # reactivate user plugins
        #     genomepy.funcions.manage_plugins("enable", active_plugins)
