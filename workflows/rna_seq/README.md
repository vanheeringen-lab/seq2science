# RNA-seq
Align samples against an assembly, count expressed genes, generate bigwigs and perform differential expression analysis. 
Bigwigs can be split for stranded protocols, and can be visualized in a trackhub or automatically generated assembly hub.  

Samples can be your own fastq files, or GSM/ERR/SRR/SRA entries. The pipeline (automatically) handles both paired-end and single-end data.
Alignment can be performed with multiple aligners. For RNA-seq, you can use STAR and Salmon (the latter requires an Ensembl assembly).
Differential gene expression (DE) analysis can optionally be performed with DESeq2.

<p align="center">
    <img src="https://raw.githubusercontent.com/vanheeringen-lab/snakemake-workflows/master/dag/rna-seq.svg?sanitize=true">
</p>

### Configuration
To specify which samples to align against which assembly the snakefile looks at the [samples.tsv](https://github.com/vanheeringen-lab/snakemake-workflows/blob/master/workflows/rna_seq/samples.tsv) file. The first entry is the column name *sample*, and every line after specifies a sample file, the second column is which assembly to align against.

In order to perform DE analysis, add columns containing the variables you with to analyse (e.g. time, treatment, batch).
Furthermore, one or more a design contrasts must be specified in the [config.yaml](https://github.com/vanheeringen-lab/snakemake-workflows/blob/master/workflows/rna_seq/config.yaml) file.

Take a look at our [rna-seq wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/5.-RNA-seq) configuration and best practices.
