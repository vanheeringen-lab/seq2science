# Workflows
## Downloading fastqs
Downloading public data in bulk from the NCBI & ENA databases has never been easier! See our [download fastq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/download_fastq) workflow.

### Pipeline steps
#### Download SRA
Most databases do not store fastq's, but they all store the raw data. For this reason for each sample submitted the raw data is downloaded. To convert this data to a fastq it has to be *dumped*. 

#### fastq-dump
The second step is translating (dumping) the raw data to a fastq file.

### Best practices
#### Filetype extensions and paired-end suffix
The pipeline will save fastq files in the *fastq_dir* directory, which is located in the *result_dir* directory by default. Fastqs already present in this directory will not be downloaded (again). 

People and tools have different preferences for storing their data. One of these differences is how to name the fastq.gz files. Some people/tools prefer `fastq.gz`, while others prefer `fq.gz`. By default the pipeline names files with the fastq.gz extension. However if you prefer to change this you can set the variable *fqsuffix* accordingly.

The same goes for paired-end suffix conventions. Some people prefer for instance `sample_pass_1`, whilst others prefer `sample_R1` (the default). You can use your preferred suffix by setting variables *fqext1* and *fqext2*.
Additionally, the pipeline will search the *fastq_dir* directory for paired-ended files with unrecognized fqexts and attempt to parse them for you. As long as the suffixes are lexicographically ordered (R1 > R2, pass1 > pass2), you can ignore the *fqext* parameters.

In the `config.yaml`:
```
fastq_dir: ./my_first_fastq_dir

fqsuffix: fq
fqext1: R1
fqext2: R2
```
or on the command line:

`snakemake --config fastq_dir=./my_first_fastq_dir fsuffix=fq fqext1=R1 fqext2=R2`


#### Downloading with ascp
ascp is a downloading protocol that allows for (sometimes much) faster speeds. If you have ascp installed you can make the workflow download through this protocol. You can specify the path to the binary in the config.yaml by key ascp_path and the ascp key in the config by key ascp_key. Or you could simply append the config on the command line, e.g.;

`snakemake --config ascp_key=$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh ascp_path=$HOME/.aspera/connect/bin/ascp`

See this [gist](https://gist.github.com/mfansler/71f09c8b6c9a95ec4e759a8ffc488be3) for easy installation of ascp.

#### Downloading many samples
Before the pipeline starts downloading it will first determine whether the samples are paired-end or single-end. This is done by doing two requests per sample with the ncbi tools. Unregistered accounts can only make 3 requests per second. This means that initialization takes +/- 0.66s per sample. If you are an extremely impatient person like me, you can up these numbers a with a [ncbi api key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) (their default one allows 10 request p/s). You can set the variables *ncbi_key* and *ncbi_requests* accordingly in the config.yaml, or on the command line, e.g.;

`snakemake --config ncbi_key=123456789 ncbi_requests=10`

The results of these lookups get cached, so when re-running parts of the workflow with the same samples no lookup online is required.

#### Fastq-dump options
We made a selection of [reasonable parameters](https://edwards.sdsu.edu/research/fastq-dump/) for the dumping of SRAs to fastq. You can change these settings by changing the variable *splot* for single-end dumping, or the variable *split* for paired-end dumping, in the config.yaml, or on the command line. 

## Alignment
Aligning samples has never been easier! See our [alignment](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/alignment) workflow.

### Pipeline steps
#### Automated trimming
The pipeline starts by trimming the reads with [trim galore!](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md). Trim galore automatically first trims the low quality 3' ends of reads, and removes short reads. After the quality trimming trim galore automatically detects which adapter was used, and trims it. The parameters of trim galore! for the pipeline can be set in the configuration by variable *trim_galore*. 

#### Alignment & Sorting
After trimming the reads are aligned against an assembly. Currently we support *bowtie2*, *bwa*, *hisat2* and *STAR* as aligners. Choosing which aligner is as easy as setting the *aligner* variable in the `config.yaml`, for example: `aligner: bwa`. Sensible defaults have been set for every aligner, but can be overwritten for either (or both) the indexing and alignment by specifying them in the `config.yaml`:
```
aligner: 
  star:
    index: '--limitGenomeGenerateRAM 55000000000`  # give STAR 55 GB RAM during indexing (default is 31 GB)
    align: '--readFilesCommand gunzip -c'          # unzips fastqs before aligning
```

The pipeline will check if the assembly you specified is present in the *genome_dir*, and otherwise will download it for you through [genomepy](https://github.com/vanheeringen-lab/genomepy). All these aligners require an index to be formed first for each assembly, but don't worry, the pipeline does this for you. 

The outputted alignment.bam is immediately sorted by either samtools or sambamba (*bam_sorter*) either in *queryname* or *coordinate* (default) order. Take a look at our [Choosing an aligner](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/3.-Alignment#choosing-an-appropriate-aligner) section for tips which aligner to use.

#### Mark duplicates
After aligning & sorting the bam duplicate reads are being 'marked' by picard Markduplicates. You can change the call by setting *markduplicates* in config.yaml.

#### Samtools index
Many downstream tools require an index of the deduplicated bam. The pipeline automatically generates these for you.

#### Quality report
Since the pipeline does a lot of steps automatically, it is very useful to have a quality report of samples. Along the way different quality control steps are taken, and are outputted in a single [multiqc report](https://multiqc.info/). Make sure to always check the report, and see our [explanation](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/9.-Interpreting-the-quality-report) on how to interpret the different scores!

#### Trackhub ###
A trackhub can be generated for this workflow. See the [Trackhub page](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/Trackhub) for more information!

### Best practices
#### (Automated) downloading of fastq
Alignment can be done with local fastq files, downloaded files, or a mix of those. Make sure to take a look at our downloading fastq [best practices](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/2.-Downloading-samples) when making use of the downloading functionality.

#### Working with local files
See our [Getting started](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/1.-Getting-started#where-does-the-workflow-store-results-and-looks-for-starting-points) about working with local files.

#### Choosing an appropriate aligner
The RNA-seq workflow supports STAR and Salmon. By default, STAR is selected.
- STAR required large amounts of memory to index the genome, but is [the more accurate](https://link.springer.com/article/10.1186/1471-2105-14-91) of the two. 
- Salmon allows for fast and light-weight transcript quantification. The accuracy can be improved by using a decoy-aware index. This can be achieved by setting *decoy_aware_index* to True in the config.yaml. However, the generation of the decoy transcripts requires an even larger amounts of memory compared to STAR.

ATAC-seq workflows often use BWA or Bowtie2. By default, BWA is selected.

#### BAM mapping quality filtering
All aligners pass a mapping quality (MAPQ) score to the reads. This score reflects the certainty that a read belongs to a certain genomic position. This score can then be used by [deeptools](https://deeptools.readthedocs.io/en/develop/) to only keep reads of which the aligner is certain where it aligns, and to limit multimappers. You can set the minimum mapq score by using the configuration variable `min_mapping_quality`. Unfortunately, most aligners use a different scoring system. This problem has been nicely summarized [here](https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/). 

For the default aligners (BWA and STAR), a minimal quality score has been set to select uniquely mapped reads only. For other aligners, we suggest the following scores to filter for uniquely mapped reads:

|Aligner|min. MAPQ scores for uniquely mapped reads|
|---|---|
|Bowtie2|41|
|BWA|30|
|Hisat2|44|
|STAR|255|

#### Cram support
For the edge-case where you want to work with cram files instead of bams, there is an option `cram_no_bam` which will convert your bam files to cram (saves around 60% storage).


## ATAC-seq
Running an ATAC-seq pipeline has never been easier! See our [ATAC-seq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/atac_seq) workflow.

### Pipeline steps
#### Peak calling

##### MACS2 
[MACS2](https://github.com/taoliu/MACS) is the successor of MACS, created by [taoliu](https://github.com/taoliu) (also one of the co-authors of hmmratac). MACS2 generates a '*pileup*' of all the reads. You can imagine this pileup as laying all reads on top of each other on the genome, and counting how high your pile gets for each basepair. MACS2 then models '*background*' values, based on the total length of all your reads and the genome length, to a [Poisson distribution](https://en.wikipedia.org/wiki/Poisson_distribution) and decides whether your experimental pileup is significant compared to the poisson distribution. MACS2 is one of the (if not the) most used peak caller.

For the calculation of peaks MACS2 requires that a genome size is being passed as one of its arguments. For the more common assemblies (e.g. mm10 and human) these numbers can be found online. However we found googling these numbers quite the hassle, so the pipeline automatically estimates this number (we calculate the genome size as **the number of unique kmers of the average read length**).  

##### Genrich
[Genrich](https://github.com/jsh58/Genrich) is the spiritual successor of MACS2, created by [John M. Gaspar](https://github.com/jsh58). Just like MACS2 is generates a '*pileup'*. However the author of genrich realized that the distribution of pileup never follows a poisson distribution. Genrich then uses a [log-normal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution) to model the background. Genrich has many nice features, but is a relatively new and untested, and is not published yet. 

##### HMMRATAC  (Not fully supported!)
[HMMRATAC](https://github.com/LiuLabUB/HMMRATAC) is a peak caller specifically made for ATAC-seq data, created by [taoliu](https://github.com/taoliu). HMMRATAC does not generate a so-called pileup. It however looks at the length distribution of reads. Since the tn5 transposase prefers to cut at nucleosome free regions in the DNA, we see that our reads are either short (at open DNA), or the length of 1 (or 2, 3, 4, etc.) nucleosomes long. HMMRATAC then classifies regions in the genome based on the length distribution of the reads in open chromatin, flanking regions, non-flanking regions and background. 

#### Replicates
All peak callers are capable of combining biological replicates. See [Replicate handling](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/4.-Replicate-handling) for more information.

#### Quality report
Since the pipeline does a lot of steps automatically, it is very useful to have a quality report of samples. Along the way different quality control steps are taken, and are outputted in a single [multiqc report](https://multiqc.info/) in the `qc` folder. Make sure to always check the report, and take a look at [interpreting the multiqc report](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/3.0-Interpreting-the-MultiQC-report)!

#### count table
A useful result the pipeline outputs is the 'count table', located at {peak_caller}/count_table_{assembly}.bed. For every peak the peak caller found in your data, it outputs the number of reads for each sample that was under that peak. 
```
				sample1		sample2
FLLO01000001.1:5649-5849	75.00000	0.00000
FLLO01000001.1:7399-7599	47.00000	1.00000
```
**Note** that this uses the number of raw reads under each peak, and does not use any kind of correction for the difference in total number of reads between samples!

#### Trackhub ###
A trackhub is be generated for this workflow. See the [Trackhub page](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/3.1-Trackhub) for more information!

### Best practices
#### (Automated) downloading of fastq
ATAC-seq can be done with local fastq files, downloaded files, or a mix of those. Make sure to take a look at our downloading fastq [best practices](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/2.-Downloading-samples) when making use of the downloading functionality.

#### Alignment
Make sure to take a look at our [alignment wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/3.-Alignment).

#### MACS2 and --shift with paired-end data
When trying to detect the cutting sites of Tn5, it is common to shift all your reads for instance 100 bp upstream, and extend them to 200 long. This makes the cutting site of a read fall in the middle of this 200 bp region. Unfortunately this is only possible in single-end mode of MACS, and makes MACS2 throw away all(!) of the paired-end mates, resulting in only half the amount of reads...

Fortunately for you, we made a script that deletes all the information that reads are paired-end after alignment from a bam. In this way we bypass the problem that shift is only possible for single-end reads. If you want to make use of this option then set the *macs2_keep_mates* variable in the configuration to True.

#### Irreproducible Discovery Rate
For idr to work properly, it needs a (large) portion of peaks that are actually not true peaks. Therefore we recommend to call peaks 'loosely'. One way of doing this is setting the [peak threshold q-value](https://github.com/taoliu/MACS#-q--pvalue) high for macs2 in the configuration.

## RNA-seq
Running an RNA-seq pipeline has never been easier! See our [RNA-seq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/rna_seq) workflow.

### Pipeline steps
#### Generate count matrices ###

In addition to aligning to the genome, [STAR](https://github.com/alexdobin/STAR) can output a reads per gene table. These counts are required for differential expression analysis. 

If alignment is no priority, either STAR or [Salmon](https://github.com/COMBINE-lab/salmon) can be used to only quantify the transcripts. This option is used when `generate_trackhub` is set to `False` (the default for this workflow).

#### Trackhub ###
A trackhub can be generated for this workflow. See the [Trackhub page](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/Trackhub) for more information!

#### Differential gene expression analysis ###
Differential expression analysis can automatically be performed using [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). Simply add one or more contrast design(s) in the `config.yaml` file. Examples are given below. Additionally, the `config.yaml` contains a commented-out set of examples. 

Note: (additional) design contrasts can be added at any time. After completing the workflow, rerunning Snakemake with new contrasts will only perform these analyses.

##### DESeq2 ####
DESeq2 automatically performs library bias correction when loading your data. Batch correction is performed if a batch effect is included in the design. After calculating differentially expressed genes, a multiple testing procedure is applied. By default this is the Benjamini-Hochberg procedure, although Independent Hypothesis Weighing can also be applied. Expression counts are log transformed (by default using the apeglm method). Finally, a list of all genes is saved to file, with analysis results for expressed genes.

***
### Contrast designs
The following section will guide you through the process of creating a DESeq2 contrast using only samples.tsv and the config.yaml. 

A contrast contains a condition, which is the variable you wish to analyse, and which contains 2 or more states. Optionally, a contrast can contain a batch effect, which you want to correct for.

#### Conditions
In samples.tsv, add a column for every property of your samples.
Examples:
* a column named 'conditions' with values ‘control’ and ‘treated’.
* a column named 'conditions' with values ‘control’, ‘treatmentA’ and ‘treatmentB’.
* a column named 'stages' with values 1, 2 and 3.

Next, in the `config.yaml` add one or more contrasts.

In order to compare two groups (e.a. control vs treated) in a column with only two groups, the contrast condition is the column name (`conditions`) (the next notation is allowed as well).

In order to compare two groups for a larger list (e.a. stage 1 vs stage 3, but not stage 2), the contrast condition is the column name, followed by the groups (`stages_3_1`).

To compare all groups against a reference (e.a. control vs both treatments A and B), the contrast condition is the column name, followed by the reference level (`condition_control`). Alternatively, group ‘all’ may also be added (`condition_all_control`).

#### Batches
If your data has batch effects, these can be accounted for. Note that this can only be done if the batch effects are spread over different conditions.

Similar to conditions, add a column to samples.tsv and note the batch of each sample:

|sample|sequencing_month|conditions|
|---|---|---|
|sample_1|jan|control|
|sample_2|feb|control|
|sample_3|jan|treatment|
|sample_4|feb|treatment|


In this case, the sequencing_month may have caused a batch effect. Since the (potential) batch effect is spread over the test conditions (control vs treatment) it can be taken into account during DE analysis. The contrast design to do this would be `sequencing_month + conditions_treatment_control` 

##### More examples

The example design `sequencing_month + conditions_treatment_control` may also be written as `sequencing_month + conditions` or `sequencing_month + conditions_all_control`. 

For consistency with regular DESeq2 contrast designs, designs may start with `~`, example: `~ sequencing_month + conditions_treatment_control`.

Our gremlins will carefully unpack all of your input and pass the unique designs to DESeq2. 

***
### Differential gene expression analysis output
Each contrast will return a table with all genes found during alignment. 

- The column `padj` contains the adjusted p-values after multiple testing. These should be used extract DE genes.
- The column `log2FoldChange` shows the fold change of each gene between the two conditions. The reference condition is the one _last mentioned_ in the contrast design, so use `conditions_treatment_control`. If you use `conditions_control_treatment` the fold change is _inverted_.
- Several other columns were kept for sake of completion, such as column `pvalue`, which contains non-adjusted p-values and _should not be used_.

Additional outputs include a blind clustering (which may reveal switched sample names), and per contrast an MA plot and a PCA plot (to accompany the differentially expressed genes table).

## ChIP-seq

Running a CHiP-seq pipeline has never been easier! See our [CHiP-seq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/chip_seq) workflow.

### Pipeline steps
#### Peak calling

##### MACS2 
[MACS2](https://github.com/taoliu/MACS) is the successor of MACS, created by [taoliu](https://github.com/taoliu) (also one of the co-authors of hmmratac). MACS2 generates a '*pileup*' of all the reads. You can imagine this pileup as laying all reads on top of each other on the genome, and counting how high your pile gets for each basepair. MACS2 then models '*background*' values, based on the total length of all your reads and the genome length, to a [Poisson distribution](https://en.wikipedia.org/wiki/Poisson_distribution) and decides whether your experimental pileup is significant compared to the poisson distribution. MACS2 is one of the (if not the) most used peak caller.

For the calculation of peaks MACS2 requires that a genome size is being passed as one of its arguments. For the more common assemblies (e.g. mm10 and human) these numbers can be found online. However we found googling these numbers quite the hassle, so the pipeline automatically estimates this number (we calculate the genome size as **the number of unique kmers of the average read length**).  

##### Genrich
[Genrich](https://github.com/jsh58/Genrich) is the spiritual successor of MACS2, created by [John M. Gaspar](https://github.com/jsh58). Just like MACS2 is generates a '*pileup'*. However the author of genrich realized that the distribution of pileup never follows a poisson distribution. Genrich then uses a [log-normal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution) to model the background. Genrich has many nice features, but is a relatively new and untested, and is not published yet. 

#### Replicates
All peak callers are capable of combining biological replicates. See [Replicate handling](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/1.1-Filling-out-the-samples.tsv#atac--and-chip-seq) for more information.

#### Quality report
Since the pipeline does a lot of steps automatically, it is very useful to have a quality report of samples. Along the way different quality control steps are taken, and are outputted in a single [multiqc report](https://multiqc.info/) in the `qc` folder. Make sure to always check the report, and take a look at [interpreting the multiqc report](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/3.0-Interpreting-the-MultiQC-report)!

#### count table
A useful result the pipeline outputs is the 'count table', located at {peak_caller}/count_table_{assembly}.bed. For every peak the peak caller found in your data, it outputs the number of reads for each sample that was under that peak. 
```
				sample1		sample2
FLLO01000001.1:5649-5849	75.00000	0.00000
FLLO01000001.1:7399-7599	47.00000	1.00000
```
**Note** that this uses the number of raw reads under each peak, and does not use any kind of correction for the difference in total number of reads between samples!

#### Trackhub ###
A trackhub is be generated for this workflow. See the [Trackhub page](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/3.1-Trackhub) for more information!

### Best practices
#### (Automated) downloading of fastq
ATAC-seq can be done with local fastq files, downloaded files, or a mix of those. Make sure to take a look at our downloading fastq [best practices](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/2.-Downloading-samples) when making use of the downloading functionality.

#### Alignment
Make sure to take a look at our [alignment wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/3.-Alignment).

#### Irreproducible Discovery Rate
For idr to work properly, it needs a (large) portion of peaks that are actually not true peaks. Therefore we recommend to call peaks 'loosely'. One way of doing this is setting the [peak threshold p-value](https://github.com/taoliu/MACS#-p--pvalue) high for macs2 in the configuration.

#### Broad peaks
Some histone marks (such as H3K9me3 and H3K27me3) produce broader (wider) peaks. If you want macs2 to call broad peaks you need to add the `--broad` flag in the macs2 call in the config. Furthermore additional tweaking of the broadness of your peaks is possible by also adding the options: `--broad-cutoff` and especially `--max-gap` in the macs call in the config.
For example using H3K27me3 data you can use  `--broad-cutoff 0.1`, and  `--max-gap 10000` (10kb), this will merge peaks closer than 10kb together as shown bellow:

WIP: Genome browser image comparing peak tracks

## scATAC-seq

Running a scATAC-seq pipeline has never been easier! See our [scATAC-seq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/scATAC_seq) workflow.

Does your scATAC protocol generate single cell fastq files? And do they need to be mapped and analysed? Then this pipeline is for your! This pipeline takes single cell fastq files, performs extensive (plate based) QC steps, and generates a binnedc SNAP object that can be used for further downstream analysis following the tutorials on: https://github.com/r3fang/SnapATAC

If instead you have a large FASTQ file containing a cell specific ID in the header, it might possible to run this pipeline with some alterations. If you have instead a Cellranger object, it is relatively easy to directly generate a SNAPobject from those, look at the doucmentation on the SnapATAC github intead.

### Pipeline steps

Concatinate samples to large FASTQ file and add cell ID to fastq header
Run Trim-Galore to remove adapters and perform QC testing of the fastq files
align using hisat2, bwa or star.

### WIP/ Optional, generate an annoation file which cell was FACsed sorted where on your plate
Depending on how you generated your scATAC data, you potentially have a FACS sorting based protocol. If this is the case you can semi-automate the annotation the barcoded fastq files to each well (and therefore FACs sorted cell type).
In order to perform this annotation, in the folder XXX WIP XXX there will be a annotation.toml file.

Do here XXX WIP XXX to generate a barcode-celltype file that can be used in downstream annotation of the SNAP-object.

### Generate and bin a SNAP_Object for downstream analysis

To find cellular heterogeneity before peak calling, SNAP-ATAC binns the genome. Using this for demensionality reduction of the data and finding cell clusters before peak calling. By default the pipeline includes bins of 5kB, but this is changable in the config.yaml file.

### How to get the pipeline started?

#### 2.6.1 Filling out the samples.tsv

#### 2.6.2 Filling out the config.yaml
