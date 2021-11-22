"""
Generate supporting files for a genome.
"""
import genomepy

genomepy.Genome(snakemake.wildcards.assembly, genomes_dir=snakemake.params.genome_dir)
