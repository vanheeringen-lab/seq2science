"""
Script to download genome annotation
"""
import os
import contextlib

import genomepy


with open(snakemake.log[0], "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):

        provider = snakemake.params.providers[snakemake.wildcards.raw_assembly]["annotation"]
        p = genomepy.ProviderBase.create(provider)

        kwargs = dict()
        if provider == "ucsc" and p.get_annotation_download_link(
                name=snakemake.wildcards.raw_assembly, kwargs={"ucsc_annotation_type": "ensembl"}):
            kwargs = {"ucsc_annotation_type": "ensembl"}

        p.download_annotation(
            name=snakemake.wildcards.raw_assembly,
            genomes_dir=snakemake.params.genome_dir,
            kwargs=kwargs)

        # sanitize the annotations
        genome = genomepy.Genome(snakemake.wildcards.raw_assembly, snakemake.params.genome_dir)
        genomepy.utils.sanitize_annotation(genome)

        # TODO remove after fixing inconsistency in gp
        if os.path.exists(snakemake.output.gtf[0][:-3]):
            genomepy.utils.gzip_and_name(snakemake.output.gtf[0][:-3])
