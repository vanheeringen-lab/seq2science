# Getting started
## Installation
TODO: conda instructions
Installation of seq2science has been made easy! All you have to do is follow these three simple steps:
- Download and install [anaconda](https://www.anaconda.com/)
  - `wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh`
  - `bash Anaconda3-2019.07-Linux-x86_64.sh` (make sure to say **yes** when asked to run conda init)
  - `source ~/.bashrc`
- clone the snakemake-workflows repository 
  - `git clone https://github.com/vanheeringen-lab/snakemake-workflows.git`
- Install the snakemake-workflows environment 
  - `conda env create -f snakemake-workflows/envs/snakemake-workflows.yaml`

## Running a workflow
TODO: update for seq2science
To run a workflow all you need to do is follow these simple steps:
- Change to a workflow directory (e.g. download_fastq)
  - `cd snakemake-workflows/workflows/download_fastq`
- Activate the (previously installed) snakemake-workflows environment 
  - `conda activate snakemake-workflows`
- Change the `config.yaml` and `samples.tsv` to your needs (see [filling out the samples.tsv](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/1.1-Filling-out-the-samples.tsv) and [filling out the config.yaml](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/1.2-filling-out-the-config.yaml) for tips.
- Run the workflow
  - `snakemake --cores 6 --use-conda`

## Where does the workflow store results and looks for 'starting points'?
Under the hood we tell snakemake which output files we expect by editing the `samples.tsv`. Snakemake then takes a look at all the rules, decides which steps it need to take to generate this output. So if for instance half of a pipeline is already run (e.g. all the samples are already aligned), snakemake will continue from this point (it will not start by trimming again). By default, Snakemake will look for all the files in the *result_dir* directory and in the appropriate folders there, and a separate directory is specified for genomes with the *genome_dir*. A typical start of a project has this layout:

```
workflow  # (e.g. alignment)
├── Snakefile
├── samples.tsv
├── config.yaml
|
├── {result_dir}
|   └── {fastq_dir}  # defaults to fastq
|       ├── sample1.fastq.gz
|       ├── sample2_pass_1.fastq.gz
|       └── sample2_pass_2.fastq.gz
|
└── {genome_dir}
    └── hg38
        └── hg38.fa
```

The organizers among us will appreciate that the directory layout can be customized to your heart's desire, as all output directory destinations can be set in the `config.yaml`. 

## Useful settings
TODO: update docs
Snakemake has a lot of settings (see [here](https://snakemake.readthedocs.io/en/stable/executable.html#all-options) for all settings), but here are a couple useful settings we recommend for starting users:
- `--use-conda`: Every rule in this repository comes with a *conda environment*. This means that for every rule a set of dependencies is specified to execute that rule successfully. When the option use-conda is used snakemake will download and install this environment for each environment. It is virtually impossible to run these workflows without this option!
- `--cores N`: The (maximum) amount of cores the workflow is allowed to use. Large parts of workflows can be executed in parallel (e.g. alignment of sample_1 and sample_2), and some applications even allow for *multithreading*. With this option you can make full use of your computer/server.
- `--keep-going`: Sometimes the workflow might encounter an error (although it really shouldn't!). With this option snakemake will continue to finish the workflow for the parts that are independent of the error.   


A typical call to the workflow then look like this:

`snakemake --cores 48 --use-conda --keep-going`

Take a look at [making a snakemake profile](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/1.3-Making-a-snakemake-profile) for making a profile in which you can specify default behaviour.


## What is Snakemake?
Since all these workflows are based on [Snakemake](https://snakemake.readthedocs.io/en/stable/)
These workflows are all based on [Snakemake](https://snakemake.readthedocs.io/en/stable/), allowing for **reproducible** and **scalable** workflows.
* Automated installation of dependencies in self-contained environments through [anaconda](https://www.anaconda.com/)
* Re-running the pipeline, either with a new configuration or new samples, is effortless
* Can be run either on personal computers or superclusters, e.g. surfsara, with minimal finetuning
* Is an extension on the Python language, which means easy to read (and maintain)



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

`snakemake file1.txt`

And snakemake will see that rule one can generate this output, sees that it requires no input, and executes the shell command. If we tell snakemake to generate file2.txt:

`snakemake file2.txt`

It will see that rule two needs to be run, takes a look at the required input for this rule, and checks whether file1.txt already exists. If it does, it will immediatly execute rule two, if file1.txt does not already exist it will execute rule one first.

Our workflows work in a way that you specify which samples you have, and which workflow you want to run, and the workflow will decide which output corresponds to the combination of those samples and workflow. Snakemake will then take a look at all the rules that are defined, and will only use the ones that are required to generate the output specific to the workflow.

For a more complete explanation of how snakemake works see the [snakemake docs](https://snakemake.readthedocs.io/).
