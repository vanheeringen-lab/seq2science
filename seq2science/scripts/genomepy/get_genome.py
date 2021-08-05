"""
Script to download genome
"""
import contextlib

import genomepy

logfile = snakemake.log[0]
providers = snakemake.params.providers
assembly = snakemake.wildcards.raw_assembly
genome_dir = snakemake.params.genome_dir
output = snakemake.output[0]

genomepy.files.rm_rf(logfile)
genomepy.logger.remove()
genomepy.logger.add(
    logfile,
    format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
    level="INFO",
)
with open(logfile, "a") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        # list user plugins
        active_plugins = genomepy.config.config.get("plugin", [])
        # deactivate user plugins
        genomepy.manage_plugins("disable", active_plugins)

        # select a provider with the annotation if possible
        a = providers[assembly]["annotation"]
        g = providers[assembly]["genome"]
        provider = g if a is None else a
        try:
            genomepy.install_genome(
                name=assembly,
                provider=provider,
                genomes_dir=genome_dir,
                force=True,
            )

            # delete the support files
            # (we recreate these separately to make the output simple for snakemake and prevent redownloading)
            genomepy.files.rm_rf(f"{output}.fai")
            genomepy.files.rm_rf(f"{output}.sizes")
            genomepy.files.rm_rf(f"{output[:-2]}gaps.bed")
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
            genomepy.manage_plugins("enable", active_plugins)
