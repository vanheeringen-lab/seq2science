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
Additionally, the pipeline will search the *fastq_dir* directory for paired-ended files with unrecognized fqexts and attempt to parse them for you. As long as the suffixes are lexicographically ordered (R1 > R2, pass1 > pass2), you can ignore the *fqext* parameters.

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
