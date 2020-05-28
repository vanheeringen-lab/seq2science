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
