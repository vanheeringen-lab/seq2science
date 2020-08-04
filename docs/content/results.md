# Using the results
## MultiQC quality report
All pipelines (except from the `download_fastq` pipeline) output a [multiQC](https://multiqc.info/) report. The report is located under `{qc_dir}/multiqc_{assembly}.html` and we highly recommend always checking out the report after a pipeline run. What is reported inside the report differs per pipeline and input.

### Rename & hide buttons
Seq2science reports quality control metrics for all your samples. However sample names on the SRA or from your own files aren't always as easy to read as you'd like. The report automatically allows for renaming of samples by the big blue buttons at the start of the report. We recommend adding a `descriptive_name` column to the samples.tsv file, where you put a human-readable sample name for each sample. This column is then automatically used by other tools as well (e.g. trackhub). 

When dealing with technical replicates and paired-end data you might find that multiqc reports way too many samples and info to properly check the report. The report also automatically generates "hide" buttons, which gives you more control over which samples are shown/hidden.   

### General statistics
The general statistics table shows a quick summary of your data. The table is interactive and you can sort it, add or remove columns (there are many hidden columns!) by clicking on the `configure columns` button.

### Metrics & Tools
* **FastQC** (raw):
  * Here the results of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are displayed **before** trimming.
* **Cutadapt**:
  * The pipeline makes use of [trimgalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for automatic adapter detection and adapter & quality trimming, which under the hood makes use of [cutadapt](https://cutadapt.readthedocs.io/en/stable/) for adapter trimming. This gives us a metric of the percentage of reads that have been trimmed.
* **FastQC** (trimmed):
  * Here the results of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are displayed **after** trimming.
* **Picard**:
  * [Picard](https://broadinstitute.github.io/picard/) is a suite of tools that can do many useful things. We use it to mark (optical and PCR) duplicates, and to get the sizes of (paired-end) inserts (useful to check for e.g. over-digestion in ATAC-seq).
* **Samtools Stats**:
  * [Samtools](http://www.htslib.org/doc/samtools-1.6.html) Stats is part of the SamTools suite. This gives us different metrics about our aligned reads.
* **deeptools**:
  * [deepTools](https://deeptools.readthedocs.io/en/develop/) is a suite of tools to process and analyze deep sequencing data. 
* **macs2 / genrich _frips**:
  * When calling peaks on your data the fraction reads in peaks score (frips) can be insightful about how well your experiment was performed. This is calculated with [featurecounts](http://subread.sourceforge.net/) of the subread module.
* **mtnucratio**:
  * In some experiments there can be a lot of contamination from mitochondrial reads and it is good to be aware of this. The `M MT/Genome reads` values give the number of reads (in millions) that are mapped to the mitochondria or genome, calculated by [MTNucRatioCalculator](https://github.com/apeltzer/MTNucRatioCalculator). You can remove mitochondrial reads with the `remove_mito` flag in the `config.yaml`.

## Trackhub
It is often good to 'eyeball' the data, and check if e.g. peak calling went alright. One of the features of the pipeline is that it can generate a trackhub, including all required support files. You can host the trackhub yourself on a web accessible location, and visualize on the [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgHubConnect). Alternatively, you can visualize the files locally in [IGV](https://software.broadinstitute.org/software/igv/).

Generation of the trackhub files is optional for all workflows, and is turned off by default. Remove `create_trackhub: False` from the config, or set the parameter to True to start generating.

### UCSC genome browser
If you move the *trackhub* folder to a web-accessible location, you can upload the URL to the `hub.txt` file on the [UCSC genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgHubConnect#unlistedHubs) to gain access to your personalized hub!

If you don't have access to a web-accessible location, the bigwig files can be manually uploaded on a UCSC trackhub, as long as the assembly used is recognized by UCSC. 

### Integrative Genomics Viewer
[IGV](https://software.broadinstitute.org/software/igv/) is a locally run genome browser with baseline functionalities for read and sequence inspection. It is an excellent alternative for quick jobs or if you do not have access to a (large enough) web-accessible location.

### BigWigs ###
Bigwigs visualize the sequencing depth per base and form the core of the trackhub. Bigwigs are stored in workflow-dependent locations, and linked in the *trackhub* folder.
- Bigwig files generated by the Alignment- and RNA-seq workflow are collected in the *bigwigs* folder. Each aligned sample (or merged sample in case of technical replicates) is converted to a bigwig file.
- Bigwig files generated by the ATAC- and ChIP-seq workflow are collected in the peak-caller directory (*macs, genrich, hmmratac*). Each condition is converted to a bigwig file.

TODO: link
See the [replicate handling page](https://github.com/vanheeringen-lab/snakemake-workflows/wiki/Replicate-handling) for more information on samples and conditions.

#### Strandedness ####
Most sequencing protocols at present are strand-specific. This specificity can be used to help identify pseudogenes originating from antisense DNA, or genes with overlapping regions on opposite strands. In order for Seq2science to use this feature, the fastq must be aligned with STAR, and the samples.tsv requires the additional column `strandedness`. This column may contain identifiers `forward` or `reverse`. If strandedness is unknown, fields may be left blank or filled with `no`.

### Genome ###
If your genome assembly is not recognized on UCSC, a number of files must be generated to map your bigwigs to. Seq2science does this for you! It creates the required *genome.2bit*, as well as a *cytobands* file. With just these files you can search your assembly by coordinates.

### Gene annotations ###
If gene annotations are available, these are added as *annotations.bigBed*. This a visible as a separate track containing genes.
Additionally, the *genome.2bit* is indexed to allow you to search your assembly by gene.

### Supporting tracks ###
Track depicting the GC-percentage and the softmasked regions of the genome are generated, similarly to [MakeHub](https://github.com/Gaius-Augustus/MakeHub).
