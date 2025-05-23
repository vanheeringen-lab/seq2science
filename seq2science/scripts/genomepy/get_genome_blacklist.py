"""
Script to download the genome blacklist.

If none exists, an empty file is generated to keep snakemake content.
"""
import os
import contextlib

import genomepy


logfile = snakemake.log[0]
fasta = snakemake.input[0]
output = snakemake.output[0]

# redirect all messages to a logfile
with open(logfile, "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        genomepy.logger.remove()
        genomepy.logger.add(
            logfile,
            format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
            level="INFO",
        )

        genome = genomepy.Genome(fasta)
        plugins = genomepy.plugins.init_plugins()
        plugins["blacklist"].after_genome_download(genome)

        # touch the file if no blacklist was available
        if os.path.exists(output):
            with open(output, "r") as f:
                lines = f.readlines()

            # only keep first three columns
            lines = ["\t".join(line.split("\t")[:3]) for line in lines]

            with open(output, "w") as f:
                f.write("\n".join(lines) + "\n")
        else:
            with open(output, 'w'):
                os.utime(output, None)
