## Differential gene/peak analysis
Seq2science outputs counts matrices for each assembly in any bulk-sequencing workflow, which can optionally can be used for differential gene/peak analysis with [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).
To do so, add one or more contrast design(s) in the `config.yaml` file, with matching identifiers in the `samples.tsv` file.

After completing the workflow, rerunning Seq2science with new contrasts will only perform the new analyses.
Alternatively, you can call the DESeq2 script from the command line.
To get started with the command line script, try out `deseq2science --help`.

This section details how the contrast designs are created.
Examples are given [below](./DESeq2.html#contrast-designs), and each `config.yaml` contains a commented-out example contrast design at the bottom.

##### single-cell DESeq2
If your counts table contains all cells, and your samples.tsv contains a group identifier for each cell, you can perform DESeq2 on this data using `deseq2science` by passing the `--single-cell` flag. 
Please note that `deseq2science` can accept a gzipped tsv file as well.

Also note that the `--single-cell` flag should *not* be used with pseudo-bulk.

##### Overview of the DESeq2 method
DESeq2 automatically performs library bias correction when loading your data, and batch correction is performed if it is included in the contrast design.
After calculating differentially expressed genes/peaks, a multiple testing procedure is applied, which is either the Benjamini-Hochberg procedure (the default) or Independent Hypothesis Weighing.
The False Discovery Rate cutoff is set by alpha, which is 0.1 by default.
Finally, count values are log transformed and shrunk (by default using the apeglm method).
These defaults can be changed in the [config.yaml](./schemas.html#deseq2), under the `deseq2` variables using the `multiple_testing_procedure`, `alpha_value` and `shrinkage_estimator` options respectively.

For more information, check out the steps in this [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), which Seq2science follows.

##### Contrast designs
The following section will guide you through the process of creating a DESeq2 contrast using only the samples.tsv and the config.yaml.

Here are some definitions: a **contrast** design tells us which samples to compare. It contains three or four parts: a **condition** (a column name in the samples.tsv), and two **groups** (labels in this column).
The first group will be the **target**, and the second group will be the **reference**.
Additionally, a contrast can optionally contain a **batch** effect, which you want to correct for.

To determine differentially expressed genes/accessible peaks, DESeq2 requires at least two samples per group.
A design contrast therefore requires at least 2x2 samples.

###### Contrast in the samples.tsv
In the `samples.tsv`, add a column for every property you wish to test in your samples.
Next, add labels to the samples involved in the test. You can leave labels empty, or add labels and not use them.
For example:

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
| sample_6 | hg38     | unused     | 3      | treatmentB |

###### Contrast in the config.yaml
Next, in the `config.yaml` add one or more contrasts:

In order to compare two groups, the contrast condition is the column name, followed by the target group and finally the reference group.
For example, to compare stage 1 to stage 3 from the examples above, the contrast would be `stages_3_1`.

To compare all groups against one reference, the contrast condition is the column name, followed by target group "all" and finally the reference group.
For example, to compare all treatments to the control from the examples above, the contrast would be `treatments_all_control`
(our gremlins will carefully unpack all of your input and pass only unique designs to DESeq2).

```
# contrasts for examples 1-3 in the config.yaml
contrasts:
  - 'conditions_knockout_wildtype'
  - 'stages_3_1'
  - 'treatments_all_control'
```

###### Group order
As a rule of thumb: the reference group should be specified **last**.

The order of the groups in a design contrast determines the direction of the expression log fold change.
In the example `stages_3_1`, a gene upregulated at stage 3 has a *positive* expression fold change.
In contrast `stages_1_3`, a gene upregulated at stage 3 has a *negative* expression fold change.
Other values are unaffected.

###### Batches
If your data has batch effects, DESeq2 can try to account for these. 
Note that this can only be done if the batch effects are spread over different conditions.

Similar to conditions, add a column to samples.tsv and note the batch of each sample:

| sample   | sequencing_month | conditions |
|----------|------------------|------------|
| sample_1 | jan              | control    |
| sample_2 | feb              | control    |
| sample_3 | jan              | treatment  |
| sample_4 | feb              | treatment  |

In this example, the `sequencing_month` may have caused a batch effect.
Since the (potential) batch effect is spread over the test conditions (control vs treatment) it can be taken into account during DE analysis.
To correct for a batch, add the column and a plus sign in from of the regular contrast.
```
contrasts:
  - 'conditions_treatment_control'
  - 'sequencing_month + conditions_treatment_control'
```

##### Output
For each contrast design, the list of *all* genes/peaks is saved to file, with analysis results for expressed genes. Briefly, these include:
- The column `padj` contains the adjusted p-values after multiple testing. **(These should be used to identify DE genes/peaks)**.
- The column `log2FoldChange` contains the fold change of each gene between the two conditions. (The reference group is the one _last mentioned_ in the contrast design, so use `condition_treatment_control`. If you use `condition_control_treatment` the fold change is _inverted_.)
- Several other columns were kept for sake of completion, such as column `pvalue`, which contains p-values not adjusted for multiple testing.

In addition, MA and PCA plots are generated for each contrast design.

If the design includes a batch effect, several PCA plots are generated to visualize the effect of the batch correction.

DESeq2 models the batch effect in their package, but downstream methods may not.
For this reason, seq2science will produce a batch-corrected counts matrix (and a batch corrected TPM matrix if quantified with Salmon).
