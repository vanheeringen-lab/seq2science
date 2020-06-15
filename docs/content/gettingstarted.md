# Getting started

## Installation

### seq2science requires anaconda

Download and install [miniconda](https://www.anaconda.com/) if not yet installed:

```console
user@comp:~$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
user@comp:~$ bash miniconda.sh # (make sure to say **yes** when asked to run conda init)
user@comp:~$ source ~/.bashrc
```

Set the correct channels (in this specific order) to use bioconda:

```console
user@comp:~$ conda config --add channels defaults
user@comp:~$ conda config --add channels bioconda
user@comp:~$ conda config --add channels conda-forge
```

### Easy installation (bioconda)

The most straightforward way to install seq2science is by using [conda](https://docs.continuum.io/anaconda/) using the [bioconda](https://bioconda.github.io/) channel. To install seq2science in a fresh environment using bioconda:

```console
(base) user@comp:~$ conda create -n seq2science seq2science
```

### Install from source

To install the latest (potentially unreleased) version of seq2science you can install from source:

```console
(base) user@comp:~$ git clone https://github.com/vanheeringen-lab/seq2science
(base) user@comp:~$ cd seq2science
(base) user@comp:~/seq2science$ conda env create --name seq2science -f requirements.yaml
(base) user@comp:~/seq2science$ conda activate seq2science
(seq2science) user@comp:~/seq2science$ pip install .
```

### Mamba

We recommend you to install [mamba](https://github.com/QuantStack/mamba) in your seq2science environment to install dependencies faster:

```console
(seq2science) user@comp:~$ conda install mamba -c conda-forge

## Running a workflow

A typical setup and run of a workflow looks like this, where you start with activating the seq2science environment.

```console
(base) user@comp:~$ conda activate seq2science
```

Then navigate to your project dir.

```console
(seq2science) user@comp:~$ cd my_project
```

Where you initialize the workflow with a configuration file and samples file, and edit those to your needs. 

```console
(seq2science) user@comp:~/my_project$ seq2science init {workflow}
```

And finally run the workflow. Note that the first time you run a workflow it will take longer as it will create the snakemake environments for all analysis steps.

```console
(seq2science) user@comp:~.my_project$ seq2science run {workflow} --cores 24
```

## Where does seq2science store results and looks for 'starting points'?

We recommend that for a typical run of seq2science you use a folder structure like this: 

```
root
└── my_project
    ├── samples.tsv
    ├── config.yaml
    └── {result_dir}
        └── {fastq_dir}
            ├── sample1.fastq.gz
            ├── sample2_R1.fastq.gz
            └── sample2_R2.fastq.gz
```

This structure is the default setting of seq2science, however can be adjusted to your liking in the config.yaml.

## What is Snakemake?

Since *under-the-hood* seq2science is based on [snakemake](https://snakemake.readthedocs.io/en/stable/), we thought it might be nice to give a little introduction to snakemake.

Snakemake is a pipeline tool that allows users to specify *rules*. Each rule is defined with what it requires as input, what it will output, and what command it needs to run to generate the output from the input. This design allows for the linking of many rules, where the input of one rule is the output of another. When invoking Snakemake it will then decide itself which rules need to be run for the output you specified. 

Here is an example Snakefile with just two (very simple) rules:

```
rule one:
    output: 
        "file1.txt"
    shell:
        # for this example we just make an empty file
        "touch {output}"

rule two:
    input:
        "file1.txt"
    output: 
        "file2.txt"
    shell:
        # for this example we just make an empty file
        "touch {output}"
```

We can tell snakemake to generate for instance file1.txt like this:

```console
(base) user@comp:~$ snakemake file1.txt
```

And snakemake will see that rule one can generate this output, sees that it requires no input, and executes the shell command. If we tell snakemake to generate file2.txt:

```console
(base) user@comp:~$ snakemake file2.txt
```

It will see that rule two needs to be run, takes a look at the required input for this rule, and checks whether file1.txt already exists. If it does, it will immediatly execute rule two, if file1.txt does not already exist it will execute rule one first.

We highly recommend everyone interested in automating a part of their analysis to take a look at snakemake! For a more complete explanation of how snakemake works see the [snakemake docs](https://snakemake.readthedocs.io/).
