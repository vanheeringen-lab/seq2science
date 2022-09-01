# Changelog
All notable changes to `seq2science` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), 
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

All changed fall under either one of these types: `Added`, `Changed`, `Deprecated`, `Removed`, `Fixed` or `Security`.

## [Unreleased]

## [0.9.5] - 2022-09-01

### Changed

- no longer writes multiqc filenames to an intermediate file
- Updated kb-python to 0.27.3

### Fixed

- downloading fastq from ena directly fixed 

## [0.9.4] - 2022-07-07

### Fixed

- (hotfix) pinned the snakemake backend for working rerunning

## [0.9.3] - 2022-06-17

### Added

- Seq2science makes a biological replicates count table, which is the mean of the biological replicates.
- Seq2science now supports differential motif analysis by gimme maelstrom!!!
- configurable setting `niceness` which sets a niceness prefix to all shell commands.

### Fixed

- issue with thread parsing when threads < 12
- seq2science should fully work with slurm
- samples moved to the cloud on SRA can be downloaded again with newest pysradb.
- issue with generating a trackhub index with tiny transcripts (we just remove them :) )
- bug with sralite files having the wrong file extension so they're not recognized as downloaded

## [0.9.2] 2022-05-30

### Added

- seq2science specific lockexception and cleanup metadata errors
- `deseq2science` now accepts the optional argument `--assembly`, which can be used if the samples.tsv contains >1 assembly to specify which one is used.
  - By default, the first assembly is used (same as before)

### Changed

- rules that download something get re-tried once, in case internet is unstable
- bam files are no longer copied when sieving is not required
- moved blacklist rules to blacklist.smk
- rule inputs now use `rules.rulename.output` where possible
- renamed `.smk` files to match the naming schemes of the other `.smk`s.
- added additional comments to clarify what happens to bam files
- cleanup cache+tarballs of conda environments, saving lots of precious disk space

### Fixed

- fixed custom assembly extensions (e.g. ERCC spike-ins) for scATAC-seq and scRNA-seq
- profiles work again
- `deseq2science` now has a clear separation between positional and optional arguments
- issue with blacklist bed containing more than 3 columns

## [0.9.1] - 2022-05-10

### Changed

- updated snakemake
  - effective genome size is now estimated per kmer length instead of per sample since checkpoints should work again.

## [0.9.0] - 2022-05-10

### Changed

- renamed most globals in uppercase (main exceptions are `config` and `samples`, `treps` and `breps`)
- moved most configuration steps into functions (reducing the number of stray globals)
- replaced static functions with dictionaries
- moved replicate stuff to the configuration
- Updated Salmon
- Added the option for Salmon to use the full genome as decoy sequence
- Salmon now uses the full genome as decoy sequence by default.
  - Config option `quantifier_decoys` controls which level of decoy aware quantification you want (options are 'none', 'partial' and 'full')
  - Option 'partial' is insanely memory intensive, and the Salmon docs suggest no benefit... 
- improved parsing of the samples.tsv. More errors early on, to prevent headache later!

### Fixed

- get_fastq_pair_reads() was using one sample, not any sample
- error message not working when trimming in scRNA-seq
- trackhubs when using a mix of stranded and unstranded datasets
- fix samples.tsv checks for forbidden symbols

## [0.8.0] - 2022-04-29

### Added

- idr call is configurable (`idr_options`)
- single-cell DESeq2 (currently only via `deseq2science` with user-specified groups per cell)
- scRNA quality control workflow with singleCellTK
  - cell calling/filtering with DropletUtils
  - mitochondrial gene set detection/filtering
  - doublet identification/filtering with scDblFinder
  - processing of alternative experiments, such as spike-in expression
  - qc report generation for cell/droplet based experiments
- added Seurat and FlatFile format export to scRNA qc workflow
- added parameter to select velocity matrix for qc and export

### Changed

- rna-seq creates a TPM table for each quantification method 
- raw/processed scRNA count tables are now stored and exported to SingleCellExperiment S4 objects instead of Seurat S4 objects 
- moved scRNA post processing to separate module
- export unspliced velocity counts to separate sce object
- seq2science should be less susceptible to poor programming environment management by using the conda-ecosystem-user-package-isolation package
- seq2science will now demand all requirements exactly the way it likes it
    - this will make the workflows more stable.
- local fastq files are no longer renamed (and should just work)
- scRNA-seq trimming code simplified

### Removed

- removed scRNA merging rule due to memory issues with large and sparse samples 
- removed deprecated scRNA post-processing workflow (superseded by singleCellTK qc workflow)

### Fixed

- fixed bug causing incorrect genome string in `read_kb_counts.R`
- bams generated with(out) filtering on size and tn5 shifting weren't removed when not necessary anymore

## [0.7.2] - 2022-03-04

### Added

- TPM to gene counts conversion with pytxi
  - by default attempts to use the GTF file to convert transcript_ids to gene_names
  - otherwise will use MyGene.info
- config option `tpm2counts` to chose which TPM to counts converter to use

### Changed

- pytxi is now the default TPM to gene counts converter (over tximeta)
- peak/gene counts tables now use descriptive names (if given)
- MultiQC DESeq2 correlation plots now display correlation metrics in the figure
- using awful practices to eliminate checkpoint strandedness
- deeptools_flags renamed to deeptools_bamcoverage
- rna-seq trackhub per base tracks by default instead of bins per 50

### Fixed

- edge-cases where seq2science was too strict with rerunning
- assembly stats log scale on the y-axis
- s2s explain wont tell you about subsampling to -1 (all) reads
- tn5 shift cigar string parsing edge-case (reference deletions/insertions)

## [0.7.1] - 2022-02-10

### Fixed

- issue with broad peaks and upsetplots

## [0.7.0] - 2022-02-02

Biggest change is that we revert back to snakemake 5.18 since higher versioned snakemake's cause instability.

### Added

- upset plot as QC for peak calling. Should give a first feeling about the distribution of peaks between samples/conditions.

### Changed

- downgraded the snakemake backend as snakemake 6+ is unstable for us.

### Fixed

- corrupt environment creation with libreadline for edgeR normalization.
- subsampling causing a crash caused by bad syntax.

## [0.6.1] - 2021-12-17

### Fixed

- corrupt environment creation with libcrypto in combination with strandedness rule

## [0.6.0] - 2021-12-12

Release 0.6.0 is a mix of bug fixes, small changes, and bigger stuff. Most importantly:

* added a deseq2science command to do differential expression analysis on user-supplied tables with seq2science settings
* for single-cell RNA-seq ADT-quantification is possible
* snakemake library updated, giving seq2science a new-ish look :) 

The full changes are listed below:

### Added

- added generic stats to the MultiQC report about the assembly, which might help pin point problems with the assembly used.
- added the slop parameter to the config.yaml of atac-seq and chip-seq workflows, just so they are more visible.
- added support for seurat object export and merging for kb workflow.
- added support for CITE-seq-count for ADT quantification
- added the option to downsample to a specific number of reads.
- new deseq2science command

### Changed

- Seq2science now makes a separate blacklist file per blacklist option (encode & mitochondria), so that e.g. RNA-seq and ATAC-seq workflows can be run in parallel and don't conflict on the blacklist.  
- error messages don't show the full traceback anymore, making it (hopefully) more clear what is going wrong.
- The effective genome size is now not calculated per sample, but per read length. When dealing with multiple samples (of similar) length this improves computational burden quite some. 
- samtools environment updated to version 1.14

### Fixed

- config option `slop` is now passed along to each rule using it
- edge-case where local samples are in the cache, but not present in the fastq_dir
- bug with differential peak/gene expression across multiple assemblies
- bug with kb ref not creating index for non-velocity analysis
- bug with count import in read_kb_counts.R for technical replicates and meta-data handling
- deseq2 ordering in multiqc report
- issue with slop not being used for the final count table
- bug with onehot peaks not reporting the sample names as columns

## [0.5.6] - 2021-10-19

### Added

- MA plot, volcano plot, and PCA plots added to QC report for deseq2 related workflows

### Changed

- updated salmon & tximeta versions
- colors for DESeq2 distance plots "fixed"
- updated bwa-mem2 version and reduced the expected memory usage of bwa-mem2 to 40GB
- seq2science now uses snakemake-minimal

### Fixed

- stranded bigwigs are no longer inverted (forward containing reverse reads and vice-versa).
- fix in `rename_sample` preventing the inversion of R1 and R2 FASTQs.
- bug with parsing cli for explanations
- show/hide buttons for treps are actually made for multiqc report
- fixes in deseq2/utils.R
  - the samples.tsv will now work with only 2 columns
  - the samples.tsv column names will be stripped of excess whitespace, similar to the config.
- ATAC-seq pipeline removing the final bams, keeping the unsorted one

## [0.5.5] - 2021-09-01

### Changed

- duplicate read marking happens before sieving and no reads get removed. Removal of duplicate reads now controlled with flag `remove_dups` in the config.
- changed option `heatmap_deeptools_options` to `deeptools_heatmap_options`
- Updated sra tools and parallel fastq-dump versions
- Updated genomepy version
- Gene annotations are no longer gzipped and ungzipped. This should reduce rerunning.

### Fixed

- rerunning being triggered too easily by input order
- issue with qc plots and broad peaks
- magic with prefetch not having the same output location on all machines
- issue with explain having duplicate lines

## [0.5.4] - 2021-07-07

### Added

- added support for kb-python kite workflow

### Changed

- kb count output validation
- optional barcodefile argument for scRNA-seq workflow
- MultiQC updated to newest version
- updated kb-python version

## [0.5.3] - 2021-06-03

### Added

- DESeq2 blind sample distance & correlation cluster heatmaps for RNA-, ATAC- ChIP-seq counts
    - find them annotated in the MultiQC when running >1 sample

### Changed

- "biological_replicate" and "technical_replicate" renamed to "..._replicates" (matches between samples.tsv & config.yaml)
- fixed bug with seq2science making a {output.allsizes} file
- Changed explain to use 'passive style'
- Genrich peak calling defaults
  - Doesn't remove PCR duplicates anymore (best to do with markduplicates)
  - Changed extsize to 200 to be similar to macs settings
  - Turned off tn5 shift, since that is done by seq2science

### Fixed

- depend less on local genomes (only when data is unavailable online)
- trackhub explanation was missing, added
- bug with broad peaks and qc that could not be made

## [0.5.2] - 2021-05-10

### Added

- added rule for scRNA post-processing R Markdown for plate/droplet based scRNA protocols (experimental)
- added explanation for kb_seurat_pp rule
- heatmap of N random peaks to the multiqc report in the end

### Fixed

- removed a warning of genome.fa.sizes already existing due to being already being downloaded beforehand (it's removed in between)
- genomepy's provider statuc checking not being used.

## [0.5.1] - 2021-04-01

### Added

- added CLI functionality to the deseq2.R script (try it with `Rscript /path/to/deseq2.R --help`!)
- --force flag to seq2science init to automatically overwrite existing samples.tsv and config.yaml
- local fastqs with Illumina's '_100' are now recognized
- added the workflow explanation to the multiqc report

### Changed

- config checks: all keys converted to lower case & duplicate keys throw an exception
- MultiQC updated to v1.10
- Link to seq2science log instead of snakemake log in final message

### Fixed

- Issue when filtering a combination of single-end and paired-end reads on template length
- explain functionality testing
- scATAC can properly use SE fastqs
- scRNA can use fqexts other than R1/R2
- fastq renaming works again
- added missing schemas to extended docs
- Bug with edgeR.upperquartile normalization. Now makes everything NaN, so pipeline finishes succesfully.

## [0.5.0] - 2021-03-03

Version 0.5.0 brings many quality of life improvements, such as seq2science automatically inferring what needs to be re-run when changing the samples.tsv and/or the config.yaml, differential peak analysis for chip/atac workflows and tab-completion!

To (hopefully) clear things up we changed the way technical and biological replicates are called, now technical and biological replicate, before replicate and condition.

It is important to note that the RNA-seq workflow DOES NOT remove duplicate reads anymore as a *default*, and that the sc/bulk ATAC-seq workflows now filters reads on the nucleosome-free region as a *default*.  

### Changed

- Keep all duplicate reads in RNA-seq by default
- Slimmed down the config printed at the start of a run
- Changed some rules into localrules when executed on a cluster
- moved onehot peaks to counts_dir
- DESeq2 contrasts now accept any column names
  - groups still cannot contain underscores
  - no longer accepts one group name
  - more examples added to the docs!

### Added

- dupRadar module to analyse read duplication types in RNA-seq
- Differential peak analysis for ATAC- and ChIP-seq!
- Options to filter bams by minimum and maximum insert sizes (added to config of bulk/sc atac)
- Support experiment ids for EBI ENA and DDBJ for downloading public samples
- More robust expression handling for BUS format detection from kb-python arguments
- Short-hand BUS syntax for indrop v1/v2
- Seq2science now supports tab-completion
- Seq2science now outputs a logfile in the directory it is run

### Fixed

- renamed more old "replicate" variables to the new "technical_replicate"
- minor logging tweak
- Chipseeker now works without defining descriptive name column
- fix bug in resources parsing of profiles
- small bug when naming a column condition in non peak-calling workflows

## [0.4.3] - 2021-01-26

### Changed

- updated tximeta to 1.6.3 and related packages to fit (now uses R 4)
- RNA-seq: sample distance matrix font scales with number of samples (should improve readability)

### Fixed

- RNA-seq: added sample distance matrix back to MulitQC
- RNA-seq: sample distance matrix legend fixed
- combine peaks with biological_replicates: keep now uses the correct peaks

## [0.4.2] - 2021-01-19

### Changed

- Updated kb-python to 0.25.1
- RNA-seq with Salmon will still use bam-related QC files if bams are generated (create_trackhub = True)

### Fixed

- gimmemotifs not working with newest pandas, now a fixed pandas version

## [0.4.1] - 2020-12-18

### Added

- more explanations for rules

### Fixed

- custom genome annotations for single-cell RNA-seq workflow
- trackhubs no longer looking for reversed strands if none are present

## [0.4.0] - 2020-12-11

### Added

- new workflow: (BETA) single cell RNA!

### Changed

- replicate renamed to technical_replicate and condition renamed to biological_replicate
- bwa-mem2 default aligner for genomic workflows, instead of bwa-mem
- interactive deeptools correlation heatmaps with static dendrograms in multiqc report
- trackhub file permissions are set to 755 so to host the files online you don't have to change those anymore

### Fixed

- bug in chip/atac trackhub generation where peaks and bigwigs used the same name, resulting in collisions and a trackhub that does not want to load 
- (literal) genome edge-case where taking the slop of peaks results in identical peaks. One of the duplicates is removed.
- IDR should work again

## [0.3.2] - 2020-11-26

### Added

- a check to see if the downloaded fastq from ENA is not empty. Related to a recent internal error (guess) at the side of ENA sending empty fastq files
- a custom message when a rule fails, that redirect to docs

### Changed

- sample layout lookup is split up in 100's, to avoid a jsondecodeerror which results from very long lists of samples
- the multiqc samples & config tables are generated in a script with its own environment to make base env smaller
- keep_mates for macs2 turned into a script with ts own environment ot make the base env smaller
- seq2science cache now respects the xdg cache
- moved genome downloading rules into scripts instead of run directives, should result in user-friendlier errors

## [0.3.1] - 2020-11-16

### Added

- Added support for multiple scrna-seq platforms (Kallistobus)
- Fastp detects the correct mate for trimming based on BUS settings.
- Support for Kallistobus short-hand syntax.

### Fixed

- make a trackhub index when the gene_name is not present in gtf file
- make a trackhub index when the gene_name is not present in gene all entries
- update Salmon & salmon rules

## [0.3.1] - 2020-11-05

### Added

- trackhub: automatic color selection
- trackhub: specify colors with the "colors" column in the samples.tsv. Accepts RGB and matplotlib colors.
- trackhub: grouped samples in a composite track with sample filters and composite control

### Changed
- updated genomepy to 0.9.1: genomes will have alternative regions removed (if designated with "alt" in the name)
- trackhub: better defaults for each track
- layouts are stored per version, as to not have collisions in the way these are stored between versions.
- scATAC no longer supports trackhub
- bigwigs are now (BPM) normalized by default

### Fixed

- markduplicates now uses $TMP_DIR, if it is defined
- RNA-seq cluster figures werent displaying text on some platforms
- not using the local annotation files
- not recognizing a mix of gzipped and unzipped annotation files
- bigwigs are now correctly labelled forward/reverse (when protocol was stranded)
- trackhub: RNA-seq trackhub now displays both strands of the bigwig (when protocol was stranded)
- trackhub: track order is now identical to the samples.tsv (was alphabetical for ChIP-/ATAC-seq)
- trackhub: assembly hub index now returns gene_name instead of transcript_id.
- bug with edgeR (upperquartile) normalization failed. Not sure why it fails, but when is does, it now returns a dataframe of nan instead of failing the rule, and thus the whole pipeline.
- use gimmemotifs 0.15.0, so gimme.combine_peaks works with numeric chromosome names
- s2s is slightly more lenient with an edge-case when running seq2science in parallel  
- clearer error message when trying samples that can not be found
- edge case with trying to dump sra from empty directory
- now give a nice error message when a technical replicate consists of a mix of paired-end and single-end samples
- issue with large number of inputs for multiqc exceeding the os command max length
- bug with downloading only SRR/DRR samples (but no GSM)
- issue with async generation of genome support files
- checking for sequencing runs when sample is already downloaded 

## [0.3.0] - 2020-09-22

### Added

- fastp as aligner (default), makes trimgalore optional other aligner
- you can now specify an url for your samples file
- RNA-seq: gene_id to gene_name conversion table will be output for downstream analysis
  - (may be empty if gtf didn't contain both fields or wrong formatting)
- RNA-seq: quantifying with salmon will now also output a gene length table
  - (gene lengths, tpms and gene counts can still be found together in the SingleCellExperiment object)

### Changed

- make use of pysradb for quering layout and SRR ids instead of API and web-scraping
- markduplicates now removes duplicates as default
- testing: clear genomepy caches between runs
- add parallel-fastq-dump fallback to fasterq-dump
- configuration rules split into more sections
- DESeq2 options renamed (from `diffexp` to `deseq2` and `contrasts`)
- DESeq2 will now generate batch corrected counts (and TPMs for Salmon) for all samples, based on the set condition column.
  - (batch corrected output is still meant for downstream analysis that cannot model batch effects independently, e.g. plotting)

### Fixed

- issue with control and technical replicates
- now also SRR numbers can be directly downloaded from ENA
- python3.8 syntaxwarnings
- chipseeker missing gtf input
- bugs with explain
- bwa-mem2 not working with less than 12 cores
- batch corrected TPMs no longer break when samples/rows are subset.

## [0.2.3] - 2020-09-01

### Changed

- retry mechanic for genomepy functions
- moved RNA-seq sample clustering to the MultiQC
- updated genomepy

### Fixed

- suffix being overwritten by layouts
- issue with combining conditions and ruleorder for macs2
- Assembly hub correctly showing annotations
- .fa.sizes staying empty

## [0.2.2] - 2020-08-24

### Added

- option to add custom files to each assembly (such as ERCC spike ins for scRNA-seq)

### Changed

- assemblies are now checked in the configuration, similar to samples
- get_genome was split in 3 rules, allowing for less reruns
- Profiles are now parsed by the s2s wrapper
- Checking for validity of samples.tsv now happens with pandasschema
- Explicit priority arguments to all group jobs (aligner + samtools_presort)
- Snakemake version (5.22.1)
- Reduced threads on salmon indexing (matching aligners)
- Make use of fasterq-dump instead of parallel-fastq-dump

### Fixed

- Test no longer use old cache files
- Profiles no longer overwrite command line arguments
- Fixed edge-case with condition column in samples but no peak-calling
- Downloading sra with prefetch tries multiple times to correct for lost connection 
- Ambiguity exception with rule narrowpeak_summit
- combine_peaks makes use of biological replicate's peaks, not technical replicate's peaks
- Bug with direct peak-calling on conditions

## [0.2.1] - 2020-08-10

### Added

- Chipseeker images in MultiQC report
- Samples that are on ENA are now directly downloaded from ENA as fastq. This means we skip the CPU instensive dumping step!

### Fixed

- Fixed issue with some samples not being findable/downloadable with s2s
- Fixed has_annotation always looking for annotation even if local files present
- Fixed bug where scatac-seq workflow was making fastqc reports per sample 

### Changed

- will try to UCSC gene annotations in Ensembl format (which uses gene IDs for the gene_id field, contrary to the UCSC format that uses transcript IDs. Wild huh?)

## [0.2.0] - 2020-08-04

### Fixed

- Allow for same condition name across different assemblies & different controls

### Added

- HISAT2 as aligner for RNA-seq
- splice-aware HISAT2 indexing for RNA-seq
- quantifier HTSeq for RNA-seq
- quantifier featurecounts for RNA-seq
- Salmon will output a gene-level TPM matrix as well
- added/expanded `seq2science explain` info (now covers RNA- and scATAC-seq too)
- sequencing strandedness may now be inferred automatically (unless specified in the config/samples.tsv)
- strandedness results are displayed in the multiQC under "Strandedness"
- a DEXSeq counts matrixs can now be generated with `dexseq: True`
- seq2science CLI now has the same reason flag as snakemake (-r/--reason flag)
- (re)added fnwi + rimls logos to the qc reports that went missing in seq2science migration

### Changed

- rules and script names in RNA-seq. ex: `txi.R` is now `quant_to_counts.R` to better reflect its function
- `quant_to_counts.R` now converts salmon transcript abundances to gene counts identically to DESeq2
- STAR no longer outputs counts, and is no longer found under `quantifiers`
- gene counts are generated from (filtered) bams when using either STAR or HISAT2 as aligner and HTSeq or featureCounts are quantifier
- batch corrected gene counts are generated if a DESeq2 design contrast inclused a batch
- batch corrected TPM are generated if a DESeq2 design contrast inclused a batch, and quantification was performed using Salmon
  - for us in ANANSE, for instance
- `seq2science explain` now retrieves messages from `explain.smk`.
- `seq2science explain` now used profiles and snakemakeOptions.

### Fixed

- the alignment workflow no longer uses strandedness
- seq2science CLI can now be run without cores with a dryrun or profile with cores
- Jenkins code style (now used mamba to install flake8)

## [0.1.0] - 2020-07-15

### Added

- bwa-mem2 as aligner
- new command-line option `explain`, which explains what has been done, and writes your material & methods section for you!

### Changed

- change the workflow names, replaced _ by -. (download_fastq to download-fastq, chip_seq to chip-seq, atac_seq to atac-seq, scatac_seq to scatac-seq, and rna_seq to rna-seq)
- changed the way seq2science is called. Moved all the logic from bin/seq2science to seq2science/cli.py

### Fixed

- Bug when merging replicates and having controls

## [0.0.3] - 2020-07-01

### Fixed
- bug when specifying 2 cores, which rounded down to zero cores for samtools sorting and crash
- edger environment was incompatible
- seq2science cache on sensible location + seq2science clean fixed
- only lookup sample layout when not local, opens up for slightly better tests in bioconda recipe

## [0.0.2] - 2020-06-29

### Fixed

- samtools using the correct nr of threads after update to v1.10

### Changed

- The count table for ATAC/ChIP-seq peaks is now made from finding all peaks within a range of 200 bp, and taking the most significant one (gimmemotifs' combine_peaks) and extending the remaining peaks 200 bp. On this count table quantile normalisation, TMM, RLE and upperquartile normalisation with CPM is done. Downstream steps log transform these and mean center them. This however means that for broadpeaks no count_table is generated.
- Snakefmt -l 121 applied

## [0.0.1] - 2020-06-17
Many minor bug- and quality of life fixes.

## [0.0.0] - 2020-06-11
First release of seq2science!

[Unreleased]: https://github.com/vanheeringen-lab/seq2science/compare/v0.9.5...develop
[0.9.5]: https://github.com/vanheeringen-lab/seq2science/compare/v0.9.4...v0.9.5
[0.9.4]: https://github.com/vanheeringen-lab/seq2science/compare/v0.9.3...v0.9.4
[0.9.3]: https://github.com/vanheeringen-lab/seq2science/compare/v0.9.2...v0.9.3
[0.9.2]: https://github.com/vanheeringen-lab/seq2science/compare/v0.9.1...v0.9.2
[0.9.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.7.2...v0.8.0
[0.7.2]: https://github.com/vanheeringen-lab/seq2science/compare/v0.7.1...v0.7.2
[0.7.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.6.1...v0.7.0
[0.6.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.5.6...v0.6.0
[0.5.6]: https://github.com/vanheeringen-lab/seq2science/compare/v0.5.5...v0.5.6
[0.5.5]: https://github.com/vanheeringen-lab/seq2science/compare/v0.5.4...v0.5.5
[0.5.4]: https://github.com/vanheeringen-lab/seq2science/compare/v0.5.3...v0.5.4
[0.5.3]: https://github.com/vanheeringen-lab/seq2science/compare/v0.5.2...v0.5.3
[0.5.2]: https://github.com/vanheeringen-lab/seq2science/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.4.3...v0.5.0
[0.4.3]: https://github.com/vanheeringen-lab/seq2science/compare/v0.4.2...v0.4.3
[0.4.2]: https://github.com/vanheeringen-lab/seq2science/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.3.2...v0.4.0
[0.3.2]: https://github.com/vanheeringen-lab/seq2science/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.2.3...v0.3.0
[0.2.3]: https://github.com/vanheeringen-lab/seq2science/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/vanheeringen-lab/seq2science/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.0.3...v0.1.0
[0.0.3]: https://github.com/vanheeringen-lab/seq2science/compare/v0.0.2...v0.0.3
[0.0.2]: https://github.com/vanheeringen-lab/seq2science/compare/v0.0.1...v0.0.2
[0.0.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.0.0...v0.0.1
[0.0.0]: https://github.com/vanheeringen-lab/seq2science/releases/tag/v0.0.0
