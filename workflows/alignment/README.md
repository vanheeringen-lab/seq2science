# Alignment
Align samples against an assembly. Samples can your own fastq files, or GSM/ERR/SRR/SRA entries. The pipeline handles (automatically) both paired-end and single-end data:

<p align="center">
    <img src="https://raw.githubusercontent.com/vanheeringen-lab/snakemake-workflows/master/dag/alignment.svg?sanitize=true">
</p>

### Configuration
To specify which samples to align against which assembly the snakefile looks at the [samples.tsv](https://github.com/vanheeringen-lab/snakemake-workflows/blob/master/workflows/alignment/samples.tsv) file. The first entry is the column name *sample*, and every line after specifies a sample file, the second column is which assembly to align against.

Take a look at our [alignment wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/Downloading-samples) configuration and best practices.
