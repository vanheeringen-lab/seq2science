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
