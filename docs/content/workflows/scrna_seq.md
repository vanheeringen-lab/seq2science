## scRNA-seq


### Pipeline steps
<p align="center">
  <img src="../../_static/scrna_seq.png">
</p>


### Downloading of sample(s)

Depending on whether the samples you start seq2science with is your own data, public data, or a mix, the pipeline might start with downloading samples. Take a look at the downloading_fastq workflow for extensive documentation about downloading of public samples.
Downloading and indexing of assembly(s)

### Downloading and indexing of assembly(s)

Depending on whether the assembly and its index you align your samples against already exist seq2science will start with downloading of the assembly through genomepy.

## Read trimming

The pipeline starts by trimming the reads with fastp!. fastp automatically detects sequencing adapters, removes short reads and performs quality filtering. The parameters of fastp! for the pipeline can be set in the configuration by variable fastp. 

**Note**: fastp currently runs in single-end mode even if you provide paired-end fastq files. This is due to the fact that reads and barcode sequences are stored in separate fastq files. After trimming, the fastq containing the reads, usually R2, may contains less reads than R1. Therefore, we perform an intermediate step by running fastq-pair to remove singleton reads before proceeding further. 



### best practices

### How to get the pipeline started?

#### 2.6.1 Filling out the samples.tsv

#### 2.6.2 Filling out the config.yaml
