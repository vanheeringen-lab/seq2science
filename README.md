# Snakemake workflows
[![Build Status](http://ocimum.science.ru.nl/jenkins/buildStatus/icon?job=Snakemake-Workflows%2Fmaster)](http://ocimum.science.ru.nl/jenkins/job/Snakemake-Workflows/job/master/)
[![star this repo](http://githubbadges.com/star.svg?user=vanheeringen-lab&repo=snakemake-workflows&style=flat)](https://github.com/vanheeringen-lab/snakemake-workflows)
[![fork this repo](http://githubbadges.com/fork.svg?user=vanheeringen-lab&repo=snakemake-workflows&style=flat)](https://github.com/boennemann/badges/fork)

The comprehensive snakemake workflows repository of the *vanheeringen lab*. Please take a look at our **[wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki)** for help with installation, how to run it, and best practices.

### Features
- Automated installation of dependencies in self-contained environments through anaconda
- Can be run either on personal computers or superclusters, e.g. surfsara, with minimal finetuning
- Automated downloading of fastq from the NCBI and ENA databases
- All our workflows work automatically with paired-end and single-end data
- Automated quality report generation
- And many more...


### Supported workflows
* [download fastq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/download_fastq) for fast and easy downloading of SRAs and parsing of fastq files from the European Nucleotide Archive and the NCBI servers.
* [alignment](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/alignment) for alignment of samples against any assembly with virtually any aligner.
* [atac seq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/atac_seq) for a complete atac-seq pipeline including peak calling and automated trackhub generation.
* [rna seq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/rna_seq) for basic gene counting and (optional) advanced differential expression analysis.
