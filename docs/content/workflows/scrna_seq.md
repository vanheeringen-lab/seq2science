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

### Filling out the config.yaml
Every workflow has many configurable options, and can be set in the config.yaml file.
In each config.yaml we highlighted a couple options that we think are relevant for that specific workflow, and set (we think) reasonable default values.

#### Custom assembly extensions
The genome and/or gene annotation can be extended with custom files, such as ERCC spike-ins for scRNA-seq.
To do so, add `custom_genome_extension: path/to/spike_in.fa` and `custom_annotation_extension: path/to/spike_in.gtf` to the config.
Seq2science will place the customized assembly in a separate folder in the `genome_dir`.
You can control the name of the customized assembly by setting `custom_assembly_suffix` in the config.

#### Quantification with kallistobus
After initializing your working directory and editing the `samples.tsv` file, specify the desired arguments for kb-pyhon via the ref (kb ref) and count (kb count) properties except for the barcode whilteist (`-w`). The path to the barcode whiltelist can be supplied via the `barcodefile` property. This step is optional since kb python python provides several pre-installed whitelists for the following technologies.

- 10XV1
- 10XV2
- 10XV3
- INDROPSV3

The white-list will be installed automatically if the appropiate technology argument is provided via the `-x` parameter in short-hand syntax.

#### BUS (Barcode/UMI/Set) format
The `-x` argument indicates the read and file positions of the UMI and barcode. Kallisto bustools should auto-detect the correct settings barcode/umi layout for the following technologies if the name is supplied:

```
name         whitelist    barcode                  umi        cDNA
---------    ---------    ---------------------    -------    -----------------------
10XV1        yes          0,0,14                   1,0,10     2,None,None
10XV2        yes          0,0,16                   0,16,26    1,None,None
10XV3        yes          0,0,16                   0,16,28    1,None,None
CELSEQ                    0,0,8                    0,8,12     1,None,None
CELSEQ2                   0,6,12                   0,0,6      1,None,None
DROPSEQ                   0,0,12                   0,12,20    1,None,None
INDROPSV1                 0,0,11 0,30,38           0,42,48    1,None,None
INDROPSV2                 1,0,11 1,30,38           1,42,48    0,None,None
INDROPSV3    yes          0,0,8 1,0,8              1,8,14     2,None,None
SCRUBSEQ                  0,0,6                    0,6,16     1,None,None
SURECELL                  0,0,6 0,21,27 0,42,48    0,51,59    1,None,None
SMARTSEQ                                                      0,None,None 1,None,None
````

Alternatively, the layout can be specified as a `bc:umi:set` triplet. The first position indicates the read, the second position the start of the feature and the third position the end of the feature. For more information and examples on the BUS format, consider the **Bus** section in the [Kallisto](https://pachterlab.github.io/kallisto/manual) manual.

###### Input preparations for KITE workflow
The steps to prepare a scRNA analysis for Feature Barcoding experiments deviates slighlty from the standard seq2science workflow. In essence, we quantify the abundance of sequence features, such as antibody barcodes, rather than a set of transcripts for a certain species. Therefore, our index does not rely on a particular assembly but is build from these sequence features. Please consider the offical [kite](https://github.com/pachterlab/kite) documentation for more details.

**1**. Prepare a two-column, tab-delimited file with your feature barcode in the first column and feature names in the second.

**Example**
|sequence|name|
|---|---|
|AACAAGACCCTTGAG|barcode 1|
|TACCCGTAATAGCGT|barcode 2|


We save this file as fb.tsv.

**2**. Copy this file to the genome folder specified in `config.yaml` where seq2science searches for assemblies.

**3**. Add the basename of the feature barcode tabel, in this case **fb**, to the assembly column in your samples.tsv file.

```
sample  assembly        
pbmc    fb      
```

An example of configuring kb-python for feature barcode analysis is shown below. Add the appropiate settings to your config.

##### Examples

Quantification (10XV3)
```
quantifier:
  kallistobus:
    count: '-x 10xv3 --h5ad --verbose'
```

RNA velocity (CEL-Seq2)
```
quantifier:
  kallistobus:
    ref: '--workflow lamanno'
    count: '-x 1,8,16:1,0,8:0,0,0 --h5ad --verbose --workflow lamanno'

barcodefile: "1col_barcode_384.tab"   
```

**Note**: The RNA velocity workflow produces count matrices for unspliced/spliced mRNA counts.  

KITE feature barcoding (CEL-Seq2)
```
quantifier:
  kallistobus:
    ref: '--workflow kite'
    count: '-x 1,8,16:1,0,8:0,0,0 --h5ad --verbose --workflow kite'

barcodefile: "1col_barcode_384.tab"    
```

#### Quantification with CITE-seq-Count
[CITE-seq-Count](https://hoohm.github.io/CITE-seq-Count/) count can be used as an alternative quantifier to pre-process ADT/Cell-hashing experiments and generate read/umi count matrices. This option cannot be used in conjunction with kallistobus.

To enable quantification with CITE-Seq-count, add the following section to your config file

Example 
```
quantifier:
  citeseqcount:
    count: '-cbf 9 -cbl 16 -umif 1 -umil 8 -cells 372 --max-error 1 --bc_collapsing_dist 1 --umi_collapsing_dist 1  -T 10 --debug'

barcodefile: "barcodes.tab"
```

#### Seurat input preparation
The seq2science scRNA workflow provides the option to automatically prepare S4 Seurat objects from kb or CITE-seq-Count workflow output. 

A Seurat object is created for each individual sample containing the raw UMI counts as default assay (RNA, ADT, spliced, unspliced). In the next step, sample-wise Seurat objects are combined and stored as a merged object. Moreover, any metadata column defined `samples.tsv` will be automatically added to each Seurat object before merging in its corresponding `@meta.data` slot. The metadata fields are spread across cell identifiers. 

All objects are stored in RDATA format and can be imported into R with the `readRDS` function. To enable Seurat object export, add the following section to your config file and adjust the Seurat object parameters depending on your analysis.

```
export_seu_objects: True

seurat_object:
    project_name: merged
    min_cells: 0
    min_features: 0
```
