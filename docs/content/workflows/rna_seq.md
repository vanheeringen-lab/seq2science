## RNA-seq
Running an RNA-seq analysis has never been easier!

### Pipeline steps
<p align="center">
  <img src="../../_static/rna_seq.png">
</p>

#### Downloading of sample(s)
Depending on whether the samples you start seq2science with is your own data, public data, or a mix, the pipeline might start with downloading samples. Take a look at the [downloading_fastq](https://vanheeringen-lab.github.io/seq2science/content/workflows/download_fastq.html) workflow for extensive documentation about downloading of public samples. 

#### Downloading and indexing of assembly(s)
Depending on whether the assembly and its index you align your samples against already exist seq2science will start with downloading of the assembly through [genomepy](https://github.com/vanheeringen-lab/genomepy).

#### Read trimming
The pipeline starts by trimming the reads with [trim galore!](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md). Trim galore automatically first trims the low quality 3' ends of reads, and removes short reads. After the quality trimming trim galore automatically detects which adapter was used, and trims it. The parameters of trim galore! for the pipeline can be set in the configuration by variable *trim_galore*. 

#### Gene counting (with HTSeq/featureCounts)
RNA-seq can be performed using gene counting methods or gene quantification methods. The following section described the former (and default) method as implemented in seq2science.

##### Alignment
Reads are aligned using `STAR` or `HISAT2` (default `STAR`). Sensible defaults have been set, but can be overwritten for either (or both) the indexing and alignment by specifying them in the `config.yaml`.

The pipeline will check if the assembly you specified is present in the *genome_dir*, and otherwise will download it for you through [genomepy](https://github.com/vanheeringen-lab/genomepy). All these aligners require an index to be formed first for each assembly, but don't worry, the pipeline does this for you.

##### Bam sieving
After aligning the bam you can choose to remove unmapped reads, low quality mappings, duplicates, and multimappers. Again, sensible defaults have been set, but can be overwritten.

##### Strandedness
Most sequencing protocols at present are strand-specific. This specificity can be used to help identify pseudogenes originating from antisense DNA, or genes with overlapping regions on opposite strands without ambiguity.
Strandedness is inferred automatically for all RNA-seq samples. For aligners it is inferred by RSeQC, the results of which can be reviewed in the MultiQC.
RSeQC inference can be overwritten by column `strandedness` in the samples.tsv. This column may contain identifiers `no`, `forward` or `reverse`. If strandedness is unknown (for some samples), fields may be left blank or filled with `nan`.
Setting `ignore_strandedness` in the config.yaml will resulting in gene counting to assume all reads are unstranded.

##### Gene counts
Gene counts are obtained from the filtered bam files using either `HTSeq` or `featureCounts` (default `HTSeq`). These counts are then combined into a count matrix per assembly for use in downstream analyses.

#### Gene quantification (with Salmon)
RNA-seq can be performed using gene counting methods or gene quantification methods. The following section described the latter method as implemented in seq2science.

#### Gene abundances
Gene abundances can be estimated using `Salmon`. Reads are aligned against the transcriptome to obtain transcript abundances (sequence strandedness is inferred automatically by Salmon), then summarized to gene-level using tximeta.
Additionally, `Salmon` generates a gene-level TPM matrix and a SingleCellExperiment object which can be opened in *R*, containing the transcript- and gene-level summaries.

***
#### Differential gene expression analysis
Seq2science outputs gene counts matrices for each assembly. Additionally, it can also perform differential expression analysis with [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).
To do so, add one or more contrast design(s) in the `config.yaml` file, with matching identifiers in the `samples.tsv` file. Examples are given below, and the RNA-seq `config.yaml` contains a commented-out set of examples.

Note: (additional) design contrasts can be added at any time. After completing the workflow, rerunning Snakemake with new contrasts will only perform these analyses.

##### DESeq2 method
DESeq2 automatically performs library bias correction when loading your data. Batch correction is performed if included in the design.
After calculating differentially expressed genes, a multiple testing procedure is applied. This is either the Benjamini-Hochberg procedure (the default) or Independent Hypothesis Weighing.
Expression counts are log transformed (by default using the apeglm method). These defaults can be changed in the `config.yaml`.

##### DESeq2 output
For each contrast design, the list of *all genes* is saved to file, with analysis results for expressed genes. Briefly, these include:
- The column `padj` contains the adjusted p-values after multiple testing. These should be used to identify DE genes.
- The column `log2FoldChange` contains the fold change of each gene between the two conditions. The reference condition is the one _last mentioned_ in the contrast design, so use `conditions_treatment_control`. If you use `conditions_control_treatment` the fold change is _inverted_.
- Several other columns were kept for sake of completion, such as column `pvalue`, which contains _non-adjusted_ p-values.

In addition, MA and PCA plots are generated for each contrast design.

If the design includes a batch effect, several PCA plots are generated to visualize the effect of the batch correction.

DESeq2 models the batch effect in their package, but downstream methods may not. For this reason, seq2science will produce a batch-corrected counts matrix (and a batch corrected TPM matrix if quantified with Salmon).

#### Differential transcript usage
Quantifying with `Salmon`, the transcript-level summaries in the SingleCellExperiment object *should* be usable for differential transcript analysis with `DEXseq`, as described [here](http://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html#dexseq).

#### Differential exon usage
Differential exon analysis by [`DEXseq`](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) can be automatically prepared by setting `dexseq: True` in the config.yaml.
This will let seq2science to output an exon counts matrix per assembly, which can be loaded directly into `DEXSeqDataSet()`.

Note: this utilizes scripts implemented by DEXseq, which are [built for Ensembl genomes](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#24_Preparing_the_annotation).

***
#### Trackhub
A UCSC compatible trackhub can be generated for this workflow. See the [trackhub page](../results.html#trackhub)<!-- @IGNORE PREVIOUS: link --> for more information!

***
### Filling out the samples.tsv
Before running a workflow you will have to specify which samples you want to run the workflow on. Each workflow starts with a `samples.tsv` as an example, and you should adapt it to your specific needs. As an example, the `samples.tsv` could look something like this:
```
sample    assembly    technical_replicate    descriptive_name    control
GSM123    GRCh38      heart_1      heart_merged        GSM234
GSM321    GRCh38      heart_1      heart_merged        GSM234
GSMabc    GRCh38      heart_2      heart_not_merged    GSM234
GSMxzy    danRer11    stage_8      stage_8             GSM234
GSM890    danRer11    stage_9      stage_9             GSM234
```

#### Sample column
This column is necessary for all workflows, not just the RNA-seq workflow. If you use the pipeline on public data this should be the name of the accession (e.g. GSM2837484). If you use the pipeline on local data this should be the *basename* of the file without the *extension(s)*. For instance, `/home/user/myfastqs/sample1.fastq.gz` would be `sample1` (for single-ended data). For paired-ended data `/home/user/myfastqs/sample2_R1.fastq.gz` and `/home/user/myfastqs/sample2_R2.fastq.gz` would be `sample2`.

#### Assembly column
This column is necessary for all workflows, except the *downloading samples* workflow. Here you simply add the name of the assembly you want your samples aligned against and the workflow will download it for you.

#### Descriptive_name column
The descriptive_name column is used for the trackhub and multiqc report. In the trackhub your tracks will be called after the descriptive name, and in the multiqc report there will be a button to rename your samples after this column. The descriptive name can not contain '-' characters, but underscores '_' are allowed.

#### Technical replicates
Technical replicates, or any fastq file you may wish to merge, are set using the `technical_replicate` column in the samples.tsv file. All samples with the same name in the `technical_replicate` column will be concatenated into one file with the replicate name.

Example `samples.tsv` utilizing replicate merging:
```
sample    assembly    technical_replicate
GSM123    GRCh38      heart
GSMabc    GRCh38      heart
GSMxzy    GRCh38      stage8
GSM890    GRCh38
```
Using this file in the alignment workflow will output *heart.bam*, *stage8.bam* and *GSM890.bam*. The MultiQC will inform you of the trimming steps performed on all samples, and subsequent information of the 'replicate' files (of which only *heart* is merged).

**Note:**
If you are working with multiple assemblies in one workflow, replicate names have to be unique between assemblies (you will receive a warning if names overlap).

#### keep
Replicate merging is turned on by default. It can be turned off by setting `technical_replicates` in the `config.yaml` to `keep`.

#### Colors column
If you are visualizing your data on the UCSC trackhub you can optionally specify the colors of each track.
To do so, you can add the color by name (google "matplotlib colors" for the options), or RGB values, in the "colors" column. Empty fields are considered black.

***
### DESeq2 contrast designs
The following section will guide you through the process of creating a DESeq2 contrast using only the samples.tsv and the config.yaml.

A contrast contains a **condition**, which is the variable you wish to analyse, and which contains 2 or more states. Optionally, a contrast can contain a **batch** effect, which you want to correct for.
To determine differentially expressed genes, DESeq2 requires at least two samples per states. A design contrast therefore requires at least 2x2 samples.

#### Conditions
In the samples.tsv, add a column for every property you wish to test in your samples. 
Next, add labels to the samples involved in the test. You can leave labels empty, or add labels and not use them. 
For example:

1. a column named 'conditions' with values ‘wildtype’ and ‘knockout’.
2. a column named 'treatments' with values ‘control’, ‘treatmentA’ and ‘treatmentB’.
3. a column named 'stages' with values 1, 2 and 3.

<table>
<tr><th>Example 1</th><th>Example 2</th><th>Example 3</th></tr>
<tr><td>

|sample|assembly|conditions|
|---|---|---|
|sample_1|hg38|wildtype|
|sample_2|hg38|knockout|
|sample_3|hg38|wildtype|
|sample_4|hg38|knockout|
|sample_5|hg38||
|sample_6|hg38|unused|

</td><td>

|sample|assembly|stages|
|---|---|---|
|sample_1|hg38|1|
|sample_2|hg38|1|
|sample_3|hg38|2|
|sample_4|hg38|2|
|sample_5|hg38|3|
|sample_6|hg38|3|

</td><td>

|sample|assembly|treatments|
|---|---|---|
|sample_1|hg38|control|
|sample_2|hg38|control|
|sample_3|hg38|treatmentA|
|sample_4|hg38|treatmentA|
|sample_5|hg38|treatmentB|
|sample_6|hg38|treatmentB|

</td></tr>
</table>

Next, in the `config.yaml` add one or more contrasts:

In order to compare two groups, the contrast condition is the column name, followed by the target group and finally the reference group.
For example, to compare stage 1 to stage 3 from the examples above, the contrast would be `stages_3_1`.

To compare all groups against one reference, the contrast condition is the column name, followed by target group "all" and finally the reference group.
For example, to compare all treatments to the control from the examples above, the contrast would be `treatments_all_control`.

```
# contrasts for examples 1-3 in the config.yaml
contrasts:
  - 'conditions_knockout_wildtype'
  - 'stages_3_1'
  - 'treatments_all_control'
```

Note: the reference group in the design determines the direction of the expression fold change.
In the example `stages_3_1`, a gene upregulated at stage 3 has a positive expression fold change.
In contrast `stages_1_3`, a gene upregulated at stage 3 has a negative expression fold change.
Other values are unaffected.

#### Batches
If your data has batch effects, DESEq2 can try to account for these. Note that this can only be done if the batch effects are spread over different conditions.

Similar to conditions, add a column to samples.tsv and note the batch of each sample:

|sample|sequencing_month|conditions|
|---|---|---|
|sample_1|jan|control|
|sample_2|feb|control|
|sample_3|jan|treatment|
|sample_4|feb|treatment|


In this example, the sequencing_month may have caused a batch effect.
Since the (potential) batch effect is spread over the test conditions (control vs treatment) it can be taken into account during DE analysis.
The contrast design to do this would be `sequencing_month + conditions_treatment_control`.

##### More design contrast examples

The example design `sequencing_month + conditions_treatment_control` may also be written as `sequencing_month + conditions_all_control`.
In each case, seq2science will understand "control" is the reference group, and compare it against the other groups in the "conditions" columns.

For consistency with regular DESeq2 contrast designs, designs may start with `~`, examples: 
`~conditions_knockout_wildtype` and `~ sequencing_month + conditions_treatment_control`.

Our gremlins will carefully unpack all of your input and pass the unique designs to DESeq2. 
