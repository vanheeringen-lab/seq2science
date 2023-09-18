## Differential gene/peak analysis
In a differential gene analysis, two groups of samples are compared in order to determine which genes are significantly over- and underexpressed in the target group, compared to the reference group.
The technique is made for RNA-seq, but can be applied to ATAC-/ChIP-seq peaks as well.

Seq2science outputs counts matrices for each assembly in any bulk-sequencing workflow, which can can be used for differential gene/peak analysis with [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).
To perform a differential analysis, add one or more contrast design(s) in the `config.yaml` file, with matching identifiers in the `samples.tsv` file.
This page details how the contrast designs are created.
Examples of design contrasts are given [below](./DESeq2.html#contrast-designs), and each `config.yaml` contains a commented-out example contrast design at the bottom.

After completing a workflow, rerunning Seq2science with new contrasts will only perform the new analyses.
Alternatively, you can call the DESeq2 script from the command line.
To get started with the command line script, try out `deseq2science --help`.

##### single-cell DESeq2
If your counts table contains all cells, and your samples.tsv contains a group identifier for each cell, you can perform DESeq2 on this data using `deseq2science` by passing the `--single-cell` flag. 
Please note that `deseq2science` can accept a gzipped tsv file as well.

Also note that the `--single-cell` flag should *not* be used with pseudo-bulk counts.

##### Overview of the DESeq2 method
DESeq2 accepts raw count data, automatically performing normalization, and batch correction if included in the contrast design.
After calculating differentially expressed genes/peaks, a multiple testing procedure is applied: either the Benjamini-Hochberg procedure (the default) or Independent Hypothesis Weighing.
The False Discovery Rate is set by the variable `alpha`, which is 0.1 by default.
Finally, count values are log transformed and shrunk (by default using the apeglm method).
These defaults can be changed in the [config.yaml](./schemas.html#deseq2), under the `deseq2` variables using the `multiple_testing_procedure`, `alpha_value` and `shrinkage_estimator` options respectively.

For more information, check out the steps in this [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), which Seq2science follows.

##### Contrast designs
The following section will guide you through the process of creating a DESeq2 contrast using only the samples.tsv and the config.yaml.

Here are some definitions: a **contrast** design tells us which groups of samples to compare. 
A contrast contains three or four parts: 

  - an optional **batch effect correction** to apply (a column in the samples.tsv)
  - a **condition** (a column in the samples.tsv)
  - two **groups** (labels in the condition column)

For each contrast design DESeq2 will determine which genes/peaks are differentially expressed in the first group (the **target** group), compared to the second group (**reference** group).
To determine differentially expressed genes/accessible peaks, DESeq2 requires at least two samples per group.
A design contrast therefore requires at least 2x2 samples.

###### Contrast in the samples.tsv
In the `samples.tsv`, add a condition column for every comparison you wish to make.
Next, add group labels to the samples involved in the test.
Here are three examples:

1. a column named 'conditions' with values 'wildtype' and 'knockout'.
2. a column named 'stages' with values 1, 2 and 3.
3. a column named 'treatments' with values 'control', 'treatmentA' and 'treatmentB'.

| sample   | assembly | conditions | stages | treatments |
|----------|----------|------------|--------|------------|
| sample_1 | hg38     | wildtype   | 1      | control    |
| sample_2 | hg38     | knockout   | 1      | control    |
| sample_3 | hg38     | wildtype   | 2      | treatmentA |
| sample_4 | hg38     | knockout   | 2      | treatmentA |
| sample_5 | hg38     |            | 3      | treatmentB |
| sample_6 | hg38     |            | 3      | treatmentB |

###### Contrast in the config.yaml
Next, in the `config.yaml` add one or more contrasts:

In order to compare two groups, the contrast condition is the column name, followed by the target group and finally the reference group.
For example: to determine which genes/peaks are differentially expressed/active at stage 3 compared to stage 1, the contrast would be `stages_3_1`.

```
contrasts:
  - stages_3_1
  - conditions_knockout_wildtype
  - treatments_all_control        # the 'all' keyword is explained below
```

The contrast `stages_3_1` will compare the expression of the target group (`3`, containing sample_5 and sample_6) against the reference group (`1`, containing sample_1 and sample_2).

The contrast `conditions_knockout_wildtype` will compare the expression of the target group (`knockout`, containing sample_2 and sample_4) against the reference group (`wildtype`, containing sample_1 and sample_3).

###### Count dispersions
For each design contrast, you can specify which samples should be used to estimate the count dispersions by giving them a label in the contrast condition column.
This can be used to exclude failed samples from the analysis.
In general, more samples improves the statistical power of the experiment.

In contrast `stages_3_1` all samples will be used to estimate the count dispersions, including sample_3 and sample_4 (because they have a label in the `stages` column).

In contrast `conditions_knockout_wildtype` sample_5 and sample_6 will **not** be used to estimate the count dispersions (because they do not have a label in the `conditions` column).

###### The 'all' keyword
To compare all groups in a column against the same reference group, the contrast condition is the column name, 
followed by the `all` keyword and finally the reference group.
For example: to compare all treatments to the control from the examples above, the contrast would be `treatments_all_control`
(our gremlins will carefully unpack all of your input and pass only unique designs to DESeq2).

###### Condition group order
As a rule of thumb: the reference group should be specified last.

The order of the groups in a design contrast determines the direction of the expression log fold change.
In the example `stages_3_1`, a gene upregulated at stage 3 has a *positive* expression fold change.
In contrast `stages_1_3`, a gene upregulated at stage 3 has a *negative* expression fold change.
Other values are unaffected.

###### Batch effect correction
If your data has batch effects, DESeq2 can try to account for these.

Similar to conditions, add a column to samples.tsv and note the batch of each sample:

| sample   | researcher | treatments |
|----------|------------|------------|
| sample_1 | jane       | control    |
| sample_2 | jane       | control    |
| sample_3 | jane       | treatmentA |
| sample_4 | jane       | treatmentA |
| sample_5 | bob        | treatmentB |
| sample_6 | jane       | treatmentB |

In this example, `jane` was sick one day and asked `bob` to perform one of her experiments. 
This may have caused a batch effect.
To correct for a batch, add the batch column, and a plus sign, in front of the regular contrast:

```
contrasts:
  - treatments_treatmentB_control                 # no batch correction
  - researcher + treatments_treatmentB_control    # batch correction
```

Note: if `bob` had performed both `treatmentB` experiments, it would be impossible to tell if the output is caused by `treatmentB` or by `bob`,
but because the batch effect is divided over the test conditions it can be taken into account during DE analysis.

##### Output
For each contrast design, the list of *all* genes/peaks is saved to file, with analysis results for expressed genes. Briefly, these include:
- The column `padj` contains the adjusted p-values after multiple testing. **(These should be used to identify DE genes/peaks)**.
- The column `log2FoldChange` contains the fold change of each gene between the two conditions.
- Several other columns were kept for the sake of completion, such as the uncorrected `pvalue`.

In addition, MA and PCA plots are generated for each contrast design (stored in the qc folder).

If a contrast design includes a batch effect, several PCA plots are generated to visualize the effect of the batch correction.
DESeq2 corrects the batch effect internally, but downstream methods may not.
Therefore, Seq2science will produce batch-corrected count and TPM matrices as well.
