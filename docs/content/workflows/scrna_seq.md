## scRNA-seq
Running a scRNA-seq pipeline has never been easier!

### Pipeline steps
<p align="center">
  <img src="../../_static/scrna_seq.png">
</p>

#### Downloading of sample(s)
Depending on whether the samples you start seq2science with is your own data, public data, or a mix, the pipeline might start with downloading samples.
You control which samples are used in the [samples.tsv](#filling-out-the-samples-tsv).
Background on public data can be found [here](./download_fastq.html#download-sra-file).

#### Downloading and indexing of assembly(s)
Depending on whether the assembly and its index you align your samples against already exist seq2science will start with downloading of the assembly through [genomepy](https://github.com/vanheeringen-lab/genomepy).

#### Read trimming
The pipeline starts by trimming the reads with fastp.
Fastp automatically detects sequencing adapters, removes short reads and performs quality filtering.
The parameters of fastp for the pipeline can be set in the configuration by the variable fastp.

**Note**: fastp currently runs in single-end mode even if you provide paired-end fastq files.
This is due to the fact that reads, barcode and umi sequences are stored in separate fastq files.
After trimming, the fastq containing the reads, usually R2, may contains less reads then R1.
Therefore, we perform an intermediate step by running `fastq-pair` to remove singleton reads before proceeding any further.

#### Quantification
Quantification is performed by running the kb-python wrapper for Kallisto bustools.
Kallisto bustools relies on pseudo-alignment of scRNA reads against a reference transcriptome index.
The resulting count matrices can be further processed with scRNA toolkits, such as Seurat or Scanpy.

### Filling out the samples.tsv
Before running a workflow you will have to specify which samples you want to run the workflow on.
Each workflow starts with a `samples.tsv` as an example, and you should adapt it to your specific needs.
As an example, the `samples.tsv` could look something like this:
```
sample  assembly        descriptive_name
pbmc    GRCh38.p13      pbmc
```

#### Sample column
If you use the pipeline on **public data** this should be the name of the accession (e.g. GSM2837484).
Accepted formats start with "GSM", "SRR", "SRX", "DRR", "DRX", "ERR" or "ERX".

If you use the pipeline on **local data** this should be the *basename* of the file without the *extension(s)*.
For example:
- `/home/user/myfastqs/sample1.fastq.gz` -------> `sample1` for single-ended data
- `/home/user/myfastqs/sample2_R1.fastq.gz` ┬> `sample2` for paired-ended data <br> `/home/user/myfastqs/sample2_R2.fastq.gz` ┘

For **local data**, some fastq files may have slightly different naming formats.
For instance, Illumina may produce a sample named `sample3_S1_L001_R1_001.fastq.gz` (and the `R2` fastq).
Seq2science will attempt to recognize these files based on the sample name `sample3`.

For **both local and public data**, identifiers used to recognize fastq files are the fastq read extensions (`R1` and `R2` by default) and the fastq suffix (`fastq` by default).
The directory where seq2science will store (or look for) fastqs is determined by the `fastq_dir` config option.
In the example above, the `fastq_dir` should be set to `/home/user/myfastqs`.
These setting can be changed in the `config.yaml`.

#### Assembly column
Here you simply add the name of the assembly you want your samples aligned against and the workflow will download it for you. In case your intend to run kb-python's kite workflow, the assembly name becomes the basename of the tab delimited feature barcode file. 

#### Descriptive_name column
The descriptive_name column is used for the multiqc report.
In the multiqc report there will be a button to rename your samples after this column.

#### Filling out the config.yaml
Every workflow has many configurable options, and can be set in the config.yaml file.
In each config.yaml we highlighted a couple options that we think are relevant for that specific workflow, and set (we think) reasonable default values.

After initializing your working directory and editing the `samples.tsv` file, specify the desired arguments for kb-pyhon via the ref (kb ref) and count (kb count) properties except for the barcode whilteist (`-w`). The path to the barcode whiltelist can be supplied via the `barcodefile` property. This step is optional since kb python python provides several pre-installed whitelists for the following technologies.

- 10XV1
- 10XV2
- 10XV3
- INDROPSV3

The white-list will be installed automatically if the appropiate technology argument is provided via the `-x` parameter in short-hand syntax.

##### BUS (Barcode/UMI/Set) format
The `-x` argument indicates the read and file positions of UMIs and barcodes in the supplied R1/R2 fastq files.
Kallisto bustools should auto-detect the correct settings if you use the short-hand syntax for your technology of choice, such as `-x 10xv2`.
Internally, this is translated to the following group of `bc:umi:set` triplets:

`0,0,16:0,16,26:1,0,0`

The ` bc:umi:set` format can be supplied as an alternative to the short-hand syntax.
For more information on the BUS format, consider the [Kallisto](https://pachterlab.github.io/kallisto/manual) manual.

##### Input preparations for KITE workflow
The steps to prepare an analysis for Feature Barcoding experiments deviate slighlty from the standard seq2science workflow. In essence, we quantify the abundance of sequence features such as antibody tags rather than transcripts. Therefore, our index does not rely on a particular assembly but is build from these sequence features. Please consider the offical [kite](https://github.com/pachterlab/kite) documentation for more details.

1. Prepare a two-column, tab-delimited (.tsv) file with your feature barcode in the first column and feature names in the second.

**Example**

|---|---|
|AACAAGACCCTTGAG|barcode 1|
|TACCCGTAATAGCGT|barcode 2|




#### Examples

##### Quantification (10XV3)
```
quantifier:
  kallistobus:
    count: '-x 10xv3 --h5ad --verbose'
```

##### RNA velocity (CEL-Seq2)
```
quantifier:
  kallistobus:
    ref: '--workflow lamanno'
    count: '-x 1,8,16:1,0,8:0,0,0 --h5ad --verbose --workflow lamanno'
```

**Note**: The RNA velocity workflow produces count matrices for unspliced/spliced mRNA.  


##### KITE feature barcoding (CEL-Seq2)
```
quantifier:
  kallistobus:
    ref: '--workflow kite'
    count: '-x 1,8,16:1,0,8:0,0,0 --h5ad --verbose --workflow kite'
```

#### Custom assembly extensions
The genome and/or gene annotation can be extended with custom files, such as ERCC spike-ins for scRNA-seq.
To do so, add `custom_genome_extension: path/to/spike_in.fa` and `custom_annotation_extension: path/to/spike_in.gtf` to the config.
Seq2science will place the customized assembly in a separate folder in the `genome_dir`.
You can control the name of the customized assembly by setting `custom_assembly_suffix` in the config.

### Filling out the config.yaml
Every workflow has many configurable options, and can be set in the `config.yaml` file.
In each `config.yaml` we highlighted a couple options that we think are relevant for that specific workflow, and set (we think) **reasonable default** values.

When a workflow starts it prints the configuration variables influencing the workflow, and (almost) all these values can be added in the `config.yaml` and changed to your liking.
You can see the complete set of configurable options in the [extensive docs](../schemas.html).

### Best practices
TODO
