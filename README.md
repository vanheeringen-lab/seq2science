# seq2science
[![Build Status](http://ocimum.science.ru.nl/jenkins/buildStatus/icon?job=Snakemake-Workflows%2Fmaster)](http://ocimum.science.ru.nl/jenkins/job/Snakemake-Workflows/job/master/lastBuild/display/redirect/)
[![star this repo](https://img.shields.io/github/stars/vanheeringen-lab/snakemake-workflows?style=flat&color=brightgreen)](https://github.com/vanheeringen-lab/snakemake-workflows/stargazers)
[![fork this repo](https://img.shields.io/github/forks/vanheeringen-lab/snakemake-workflows?style=flat&color=brightgreen)](https://github.com/vanheeringen-lab/snakemake-workflows/fork)

The comprehensive snakemake workflows repository of the *van heeringen lab*. Please take a look at our [docs](https://github.com/vanheeringen-lab/snakemake-workflows/wiki) for help with installation, how to run it, and best practices.

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
* [chip seq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/chip_seq) for a complete CHiP-seq pipeline including peak calling and automated trackhub generation.
