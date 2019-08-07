# Download fastq
Download fastq files from the ENA or NCBI servers. 

### Configuration
To specify which samples to download the snakefile looks at the [samples.tsv](https://github.com/vanheeringen-lab/snakemake-workflows/blob/master/workflows/download_fastq/samples.tsv) file. The first entry is the column name *sample*, and every line after specifies a sample to be downloaded. Samples are downloaded from the NCBI and ENA databases. 

Take a look at our [downloading samples wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/Downloading-samples) for best practices.
