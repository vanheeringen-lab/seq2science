# Filling out the samples.tsv
Before running a workflow you will have to specify which samples you want to run the workflow on. Each workflow starts with a `samples.tsv` as an example, and you should adapt it to your specific needs. As an example `samples.tsv` could look like this, but not all columns are always necessary or required:

```
sample    assembly    replicate    descriptive_name
GSM123    GRCh38      heart_1      heart_merged
GSM321    GRCh38      heart_1      heart_merged
GSMabc    GRCh38      heart_2      heart_not_merged
GSMxzy    danRer11    stage_8      stage_8
GSM890    danRer11    stage_9      stage_9
```

## sample column
This column is necessary for all workflows. If you use the pipeline on public data this should be the name of the accession (e.g. GSM2837484), if you use the pipeline on local data this should be the *basename* of the file without the *extension*. For instance `/home/user/mfastqs/sample1.fastq.gz` would be sample1.

## assembly column
This column is necessary for all workflows, except the *downloading samples* workflow. Here you simply add the name of the assembly you want your samples aligned against and the workflow will download it for you.

## descriptive_name column
The descriptive_name column is used for the trackhub and multiqc report. In the trackhub your tracks will be called after the descriptive name, and in the multiqc report there will be a button to rename your samples after this column.

## working with technical and biological replicates
### Technical replicates
Technical replicates, or any fastq file you may wish to merge, are set using the `replicate` column in the samples.tsv file. All samples with the same name in the `replicate` column will be concatenated into one file with the replicate name.

Example samples.tsv utilizing replicate merging:
```
sample    assembly    replicate
GSM123    GRCh38      heart
GSMabc    GRCh38      heart
GSMxzy    danRer11    stage8
GSM890    danRer11
```
Using this file in the alignment workflow will output *heart.bam*, *stage8.bam* and *GSM890.bam*. The MultiQC will inform you of the trimming steps performed on all samples, and subsequent information of the 'replicate' files (of which only *heart* is merged).

#### Notes ####
If you are working with multiple assemblies in one workflow, replicate names have to be unique between assemblies (you will receive a warning if names overlap).

Replicate merging is turned on by default. It can be turned off by setting `technical_replicates` to `keep`.

### Biological replicates
#### RNA-seq
During the RNA-seq workflow, biological replicates are specified in the contrast design, as described in the [RNA-seq chapter](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/6.-RNA-seq).

#### ATAC- and ChIP-seq
During ATAC- and ChIP-seq workflows, peak calling can be performed per biological condition, depending on your configuration setting. Biological conditions are determined by the `condition` column in the samples.tsv file. How these samples are handled is specified by configuration variable `biological_replicates`.

If you merge technical replicates *and* combine biological replicates, the technical replicates are merged first.

#### Keep
The default setting. All samples/merged replicates are analyzed individually. The `condition` column is ignored.

#### Irreproducible Discovery Rate (IDR)
One of the more common methods to combine biological replicates is by the irreproducible discovery rate ([idr](https://github.com/kundajelab/idr)). Shortly; idr sorts all the peaks of two replicates separately on their significance. Since true peaks should be very significant for both replicates these peaks will be one of the highest sorted peaks. As peaks get less and less true (and thus their significance) their ordering also becomes more random between the samples. The idr method then only keeps the peak that overlap *nonrandomly* between the samples. The idr method only works for two replicates, so can not be used when you have more than 2 (`n == 2`).

#### Fisher's method
[Fisher's method ](https://en.wikipedia.org/wiki/Fisher%27s_method) simply is a method to 'combine' multiple p-values. This method is built-in for genrich and macs2, and allows any number of replicates (`n >= 2`). However a 'disadvantage' of this method is that it assumes the p-values a method generates are actual p-values (and not a general indication of significance). MACS2 is slightly notorious for that its p-values are not true p-values. Genrich claims by fitting on a log-normal distribution that their p-values are true. Even though fisher's method might not be ideal, it's the only option you have (we provide) within this pipeline for more than 2 replicates.

## Final notes
- Make sure that the samples.tsv is a tab separated values file when running the pipeline.
- Feel free to delete or add columns to your liking.