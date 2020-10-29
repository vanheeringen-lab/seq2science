"""
add description
"""
import genomepy

genomepy.Genome(snakemake.wildcards.assembly, genomes_dir=config["genome_dir"])
