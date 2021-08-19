"""
Script to download genome annotation
"""
import contextlib

import genomepy


with open(snakemake.log[0], "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        # list user plugins
        active_plugins = genomepy.functions.config.get("plugin", [])
        # deactivate user plugins
        genomepy.functions.manage_plugins("disable", active_plugins)

        # select a provider with the annotation if possible
        provider = snakemake.params.providers[snakemake.wildcards.raw_assembly]["annotation"]

        kwargs = dict()
        if provider == "ucsc":
            p = genomepy.ProviderBase.create(provider)
            kwargs = {"ucsc_annotation_type": "ensembl"}
            if not p.get_annotation_download_link(
                    name=snakemake.wildcards.raw_assembly,
                    kwargs=kwargs
            ):
                kwargs = dict()

        try:
            genomepy.install_genome(
                name=snakemake.wildcards.raw_assembly,
                provider=provider,
                genomes_dir=snakemake.params.genome_dir,
                only_annotation=True,
                force=True,
                kwargs=kwargs,
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
            # reactivate user plugins
            genomepy.functions.manage_plugins("enable", active_plugins)
