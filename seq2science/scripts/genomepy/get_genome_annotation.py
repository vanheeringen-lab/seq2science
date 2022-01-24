"""
Script to download genome annotation
"""
import contextlib

import genomepy


logfile = snakemake.log[0]
assembly = snakemake.wildcards.raw_assembly
providers = snakemake.params.providers
provider = snakemake.params.provider
genome_dir = snakemake.params.genome_dir

# redirect all messages to a logfile
with open(logfile, "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        genomepy.logger.remove()
        genomepy.logger.add(
            logfile,
            format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
            level="INFO",
        )

        # list user plugins
        active_plugins = genomepy.config.config.get("plugin", [])
        if active_plugins:
            print("Deactivating user plugins")
            genomepy.manage_plugins("disable", active_plugins)

        # select user specified provider
        if provider is None:
            # select a provider with the annotation
            provider = providers[assembly]["annotation"]

        try:
            genomepy.install_genome(
                name=assembly,
                provider=provider,
                genomes_dir=genome_dir,
                only_annotation=True,
                force=True,
            )

        except Exception as e:
            print(e)
            print("\nSomething went wrong while downloading the gene annotation (see error message above). "
                  "When this happens it is almost always because we had troubles connecting to the"
                  "servers hosting the genome assemblies. Usually this is resolved by just running seq2science"
                  "again, either immediately or in a couple hours.\n\n"
                  "If the problem persists you could try running `seq2science clean` and see if that resolves the "
                  "issue.")

        finally:
            if active_plugins:
                print("Reactivating user plugins")
                genomepy.manage_plugins("enable", active_plugins)
