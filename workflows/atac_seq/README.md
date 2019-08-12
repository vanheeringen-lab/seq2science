# ATAC-seq
Call peaks and generate a trackhub. Samples can your own fastq files, or GSM/ERR/SRR/SRA entries. The pipeline (automatically) handles both paired-end and single-end data.

TODO figure

### Configuration
To specify which samples to run through the pipeline the snakefile looks at the [samples.tsv](https://github.com/vanheeringen-lab/snakemake-workflows/blob/master/workflows/atac_seq/samples.tsv) file. The first entry is the column name *sample*, and every line after specifies a sample to be downloaded. The second column is which *assembly* to align and peak call against. An optional third column exists which specifies the experimental condition which allows the pipeline to combine replicates through IDR.

Take a look at our [ATAC-seq wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/4.-ATAC-seq) for best practices.
