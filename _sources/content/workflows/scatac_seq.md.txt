## scATAC-seq
Running a scATAC-seq pipeline has never been easier! See our [scATAC-seq](https://github.com/vanheeringen-lab/snakemake-workflows/tree/master/workflows/scATAC_seq) workflow.

Does your scATAC protocol generate single cell fastq files? And do they need to be mapped and analysed? Then this pipeline is for your! This pipeline takes single cell fastq files, performs extensive (plate based) QC steps, and generates a binned  SNAP-object that can be used for further downstream analysis following the tutorials on: https://github.com/r3fang/SnapATAC

If instead you have a large FASTQ file containing a cell specific ID in the header, it might possible to run this pipeline with some alterations. If you have instead a Cellranger object, it is relatively easy to directly generate a SNAPobject, look at the documentation on the SnapATAC github intead.

### Pipeline steps
#### WIP: (Automated) downloading of fastq
scATAC-seq can be done with local fastq files, downloaded files, or a mix of those. Make sure to take a look at our downloading fastq [best practices](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/2.-Downloading-samples) when making use of the downloading functionality.

#### merge cell fastq files based on processing  
Single cell fastq files get a barcode id corresponding to the sample name. This is added to the fastq header, after which each single cell fastq is merged to a large FASTQ file with all the technical replica's (for example all cells from a plate). This results in a multi-QC report containing QC per technical replica's. This makes lower plate quality easy to spot.
Additional single cell QC will be available in downstream analysis of the SNAP-object.
See [Replicate handling](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/4.-Replicate-handling) for more information.
See [Filling out the samples.tsv and config.yaml](https://github.com/vanheeringen-lab/snakemake-workflows/blob/docs/docs/fillingout.md) for more information.

#### Alignment
Make sure to take a look at our [alignment wiki](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/3.-Alignment).

#### Generate and bin a SNAP_Object for downstream analysis
To find cellular heterogeneity before peak calling, SNAP-ATAC binns the genome. Using this for demensionality reduction of the data and finding cell clusters before peak calling. By default the pipeline includes bins of 5kB, but this is changable in the config.yaml file.

### best practices
TODO: Make sure to take a look at: [example preprocessing](scATAC_postprocessing.html)

### How to get the pipeline started?

#### 2.6.1 Filling out the samples.tsv

#### 2.6.2 Filling out the config.yaml
