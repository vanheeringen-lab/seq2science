# Download fastq
Download fastq files from the ENA or NCBI servers as SRA and (parallel) dump them as fastq. The pipeline can handle both paired-end and single-end data:

<p align="center">
    <img src="https://raw.githubusercontent.com/vanheeringen-lab/snakemake-workflows/master/dag/download_fastq.svg?sanitize=true">
</p>

### Configuration
To specify which samples to download the snakefile looks at the [samples.tsv](https://github.com/vanheeringen-lab/snakemake-workflows/blob/master/workflows/download_fastq/samples.tsv) file. The first entry is the column name *sample*, and every line after specifies a sample to be downloaded. Samples are downloaded from the NCBI and ENA databases. 

Take a look at our [downloading samples wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/2.-Downloading-samples) for best practices.
