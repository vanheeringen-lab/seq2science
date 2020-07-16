## RNA-seq
Running an RNA-seq analysis has never been easier!

### Pipeline steps
<p align="center">
  <img src="../../_static/rna_seq.png">
</p>

#### Downloading of sample(s)
Depending on whether the samples you start seq2science with is your own data, public data, or a mix, the pipeline might start with downloading samples. Take a look at the [downloading_fastq](download_fastq.md) workflow for extensive documentation about downloading of public samples. 

#### Downloading and indexing of assembly(s)
Depending on whether the assembly and its index you align your samples against already exist seq2science will start with downloading of the assembly through [genomepy](https://github.com/vanheeringen-lab/genomepy).

#### Read trimming
The pipeline starts by trimming the reads with [trim galore!](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md). Trim galore automatically first trims the low quality 3' ends of reads, and removes short reads. After the quality trimming trim galore automatically detects which adapter was used, and trims it. The parameters of trim galore! for the pipeline can be set in the configuration by variable *trim_galore*. 

#### Quantification & gene count matrices
Gene counts can be obtained by two distinct methods: either by directly quantifying transcript abundances using `Salmon`, or by summarizing counts from bam files.
For the latter approach, fastqs are aligned by splice-aware aligners `STAR` or `HISAT2`. Next, the bam files are filtered according to configurable specification (which does not happen for Salmon), and counts are quantified by either `HTSeq-count` or `featureCounts`.

Gene counts are aggregated per assembly into a count matrix. Additionally, `salmon` generates a SingleCellExperiment object which can be opened in *R*, containing the transcript- and gene-level summaries. 

***
### Optional: Trackhub
A UCSC compatible trackhub can be generated for this workflow. See the [trackhub page](../results.html#trackhub) for more information!

#### Alignment
To generate a trackhub the trimmed reads are aligned against an assembly using `star` or `HISAT2`. Sensible defaults have been set, but can be overwritten for either (or both) the indexing and alignment by specifying them in the `config.yaml`.

The pipeline will check if the assembly you specified is present in the *genome_dir*, and otherwise will download it for you through [genomepy](https://github.com/vanheeringen-lab/genomepy). All these aligners require an index to be formed first for each assembly, but don't worry, the pipeline does this for you. 

#### Bam sieving
After aligning the bam you can choose to remove unmapped reads, low quality mappings, duplicates, and multimappers.

#### Strandedness
Most sequencing protocols at present are strand-specific. This specificity can be used to help identify pseudogenes originating from antisense DNA, or genes with overlapping regions on opposite strands. In order for seq2science to use this feature, the fastq must be aligned with STAR, and the samples.tsv requires the additional column `strandedness`. This column may contain identifiers `forward` or `reverse`. If strandedness is unknown, fields may be left blank or filled with `no`.

***
### Optional: Differential gene expression analysis
Differential expression analysis can automatically be performed using [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). Simply add one or more contrast design(s) in the `config.yaml` file. Examples are given below. Additionally, the `config.yaml` contains a commented-out set of examples. 

Note: (additional) design contrasts can be added at any time. After completing the workflow, rerunning Snakemake with new contrasts will only perform these analyses.

##### DESeq2
DESeq2 automatically performs library bias correction when loading your data. Batch correction is performed if included in the design. 
After calculating differentially expressed genes, a multiple testing procedure is applied. This is either the Benjamini-Hochberg procedure (the default) or Independent Hypothesis Weighing. 
Expression counts are log transformed (by default using the apeglm method). These defaults can be changed in the `config.yaml`. 
Finally, the list of all genes is saved to file, with analysis results for expressed genes.

In addition, MA and PCA plots are generated for each contrast design. If the design includes a batch effect, several PCA plots are generated to visualize the effect of the batch correction.

DESeq2 models the batch effect in their package, but downstream methods may not. For this reason, seq2science will produce a batch-corrected counts matrix (and a batch corrected TPM matrix if quantified with Salmon). 

***
### Filling out the samples.tsv
Before running a workflow you will have to specify which samples you want to run the workflow on. Each workflow starts with a `samples.tsv` as an example, and you should adapt it to your specific needs. As an example, the `samples.tsv` could look something like this:
```
sample    assembly    replicate    descriptive_name    control
GSM123    GRCh38      heart_1      heart_merged        GSM234
GSM321    GRCh38      heart_1      heart_merged        GSM234
GSMabc    GRCh38      heart_2      heart_not_merged    GSM234
GSMxzy    danRer11    stage_8      stage_8             GSM234
GSM890    danRer11    stage_9      stage_9             GSM234
```

#### Sample column
This column is necessary for all workflows, not just the atac-seq workflow. If you use the pipeline on public data this should be the name of the accession (e.g. GSM2837484), if you use the pipeline on local data this should be the *basename* of the file without the *extension*. For instance `/home/user/myfastqs/sample1.fastq.gz` would be `sample1`.

#### Assembly column
This column is necessary for all workflows, except the *downloading samples* workflow. Here you simply add the name of the assembly you want your samples aligned against and the workflow will download it for you.

#### Descriptive_name column
The descriptive_name column is used for the trackhub and multiqc report. In the trackhub your tracks will be called after the descriptive name, and in the multiqc report there will be a button to rename your samples after this column.

#### Technical replicates
Technical replicates, or any fastq file you may wish to merge, are set using the `replicate` column in the samples.tsv file. All samples with the same name in the `replicate` column will be concatenated into one file with the replicate name.

Example `samples.tsv` utilizing replicate merging:
```
sample    assembly    replicate
GSM123    GRCh38      heart
GSMabc    GRCh38      heart
GSMxzy    danRer11    stage8
GSM890    danRer11
```
Using this file in the alignment workflow will output *heart.bam*, *stage8.bam* and *GSM890.bam*. The MultiQC will inform you of the trimming steps performed on all samples, and subsequent information of the 'replicate' files (of which only *heart* is merged).

**Note:**
If you are working with multiple assemblies in one workflow, replicate names have to be unique between assemblies (you will receive a warning if names overlap).

Replicate merging is turned on by default. It can be turned off by setting `technical_replicates` in the `config.yaml` to `keep`.

#### Strandedness ####
To split bigwigs by strand, the samples.tsv requires the column `strandedness`. This column may contain identifiers `forward` or `reverse`. If strandedness is unknown, fields may be left blank or filled with `no`.


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

In order to compare two groups from a list (e.a. stage 1 vs stage 3), the contrast condition is the column name, followed by the groups (`stages_3_1`). Note that the second group (stage 1 in this example) is the reference group. 

To compare all groups against a reference (e.a. control vs both treatments A and B), the contrast condition is the column name, followed by the reference group (`condition_control`). Alternatively, group ‘all’ may also be added (`condition_all_control`).

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

##### More design contrast examples

The example design `sequencing_month + conditions_treatment_control` may also be written as `sequencing_month + conditions_all_control`, `sequencing_month + conditions_control_all` or `sequencing_month + conditions_control`. In each case, seq2science will understand "control" is the reference group, and compare it agains the other groups in the "conditions" columns. 

For consistency with regular DESeq2 contrast designs, designs may start with `~`, example: `~ sequencing_month + conditions_treatment_control`.

Our gremlins will carefully unpack all of your input and pass the unique designs to DESeq2. 

***
### Differential gene expression analysis output
Each contrast will return a table with all genes found during alignment. 

- The column `padj` contains the adjusted p-values after multiple testing. These should be used extract DE genes.
- The column `log2FoldChange` shows the fold change of each gene between the two conditions. The reference condition is the one _last mentioned_ in the contrast design, so use `conditions_treatment_control`. If you use `conditions_control_treatment` the fold change is _inverted_.
- Several other columns were kept for sake of completion, such as column `pvalue`, which contains non-adjusted p-values.

Additional outputs include a blind clustering (which may reveal switched sample names), and per contrast an MA plot and a PCA plot (to accompany the differentially expressed genes table).
