## scRNA-seq
Running a scRNA-seq pipeline has never been easier!

### Pipeline steps
<p align="center">
  <img src="../../_static/scrna_seq.png">
</p>

#### Downloading of sample(s)
Depending on whether the samples you start seq2science with is your own data, public data, or a mix, the pipeline might start with downloading samples.
Take a look at the [downloading_fastq](./download_fastq.html) workflow for extensive documentation about downloading of public samples.

### Downloading and indexing of assembly(s)

Depending on whether the assembly and its index you align your samples against already exist seq2science will start with downloading of the assembly through genomepy.

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

### Best practices

### How to get the pipeline started?

### Filling out the samples.tsv

Before running a workflow you will have to specify which samples you want to run the workflow on.
Each workflow starts with a `samples.tsv` as an example, and you should adapt it to your specific needs.
As an example, the `samples.tsv` could look something like this:

```
sample  assembly        descriptive_name
pbmc    GRCh38.p13      pbmc
```

#### Sample column
If you use the pipeline on public data this should be the name of the accession (e.g. GSM2837484).
(Accepted formats start with "GSM", "SRR", "SRX", "DRR", "DRX", "ERR" or "ERX")

If you use the pipeline on local data this should be the *basename* of the file without the *extension(s)*.
For instance, `/home/user/myfastqs/sample1.fastq.gz` would be `sample1` (for single-ended data).
For paired-ended data `/home/user/myfastqs/sample2_R1.fastq.gz` and `/home/user/myfastqs/sample2_R2.fastq.gz` would be `sample2`.

Some local fastq files may have slightly different naming formats.
Illumina's sample for instance may produce a paired-ended sample with these fastqs: `sample3_S1_L001_R1_001.fastq.gz` and `sample3_S1_L001_R2_001.fastq.gz`.
Seq2science will attempt to recognize these files based on the sample name (in this case `sample3`).
Other identifiers used are the fastq read extensions (`R1` and `R2` by default) and the fastq suffix (`fastq` by default), which can be changed in the config.

#### Assembly column
Here you simply add the name of the assembly you want your samples aligned against and the workflow will download it for you.

#### Descriptive_name column
The descriptive_name column is used for the multiqc report.
In the multiqc report there will be a button to rename your samples after this column.

#### Filling out the config.yaml

Every workflow has many configurable options, and can be set in the config.yaml file.
In each config.yaml we highlighted a couple options that we think are relevant for that specific workflow, and set (we think) reasonable default values.

After initializing your working directory and editing the `samples.tsv` file, you have to specify if you either want to perform quantification or RNA velocity analysis.
For velocity analysis, add the `--workflow lamanno` argument to the ref and count properties as shown below.

##### Quantifier

```
quantifier:
  kallistobus:
    count: '-x 10xv3 --h5ad --verbose'
```

##### RNA velocity 
```
quantifier:
  kallistobus:
    ref: '--workflow lamanno'
    count: '-x 10xv3 --h5ad --verbose --workflow lamanno'
```

This will ensure that kallisto bustools generates the required files for RNA velocity estimation and produces count matrices for unspliced/spliced mRNA.  

##### BUS (Barcode/UMI/Set) format

The `-x` argument indicates the read and file positions of UMIs and barcodes in the supplied R1/R2 fastq files.
Kallisto bustools should auto-detect the correct settings if you use the short-hand syntax for your technology of choice, such as `-x 10xv2`.
Internally, this is translated to the following group of `bc:umi:set` triplets:

`0,0,16:0,16,26:1,0,0`

The ` bc:umi:set` format can be supplied as an alternative to the short-hand syntax.
For more information on the BUS format, consider the [Kallisto](https://pachterlab.github.io/kallisto/manual) manual.

#### Custom assembly extensions
The genome and/or gene annotation can be extended with custom files, such as ERCC spike-ins for scRNA-seq.
To do so, add `custom_genome_extension: path/to/spike_in.fa` and `custom_annotation_extension: path/to/spike_in.gtf` to the config.
Seq2science will place the customized assembly in a separate folder in the `genome_dir`.
You can control the name of the customized assembly by setting `custom_assembly_suffix` in the config.
