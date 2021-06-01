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
