# Welcome
The time of hacking shell scripts into a pipeline, less stable than the average house of cards, is no more! This repository is the attempt of the *van heeringen lab* to generate an all encompassing collection of pipelines which can be used by complete beginners to bioinformatics and experienced bioinformaticians alike. 

Automate the boring stuff with ~~Python~~ [Snakemake](https://snakemake.readthedocs.io/en/stable/)! These workflows are all based on snakemake, allowing for **reproducible** and **scalable** workflows. 
* Automated installation of dependencies in self-contained environments through [anaconda](https://www.anaconda.com/)
* Re-running the pipeline, either with a new configuration or new samples, is effortless
* Can be run either on personal computers or superclusters, e.g. surfsara, with minimal finetuning
* Is an extension on the Python language, which means easy to read (and maintain)

Some of seq2science's cool features:
* Automated **downloading of fastq** from the NCBI and ENA databases.
* Automated **downloading of genomes** through genomepy.
* All our workflows work automatically with **paired-end** and **single-end** data.
* Most workflows allow for a **wide range of tools** to be used.
* Automated **quality report** generation.
* Automated **trackhub** generation.
* And many more...

## First time users
Please take a look at our [getting started](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/1.-Getting-started) section.
