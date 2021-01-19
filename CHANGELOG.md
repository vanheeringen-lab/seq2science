# Changelog
All notable changes to `seq2science` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), 
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

All changed fall under either one of these types: `Added`, `Changed`, `Deprecated`, `Removed`, `Fixed` or `Security`.

## [Unreleased]

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

[Unreleased]: https://github.com/vanheeringen-lab/seq2science/compare/master...v0.4.1
[0.4.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.4.1...v0.4.0
[0.4.0]: https://github.com/vanheeringen-lab/seq2science/compare/v0.4.0...v0.3.2
[0.3.2]: https://github.com/vanheeringen-lab/seq2science/compare/v0.3.2...v0.3.1
[0.3.1]: https://github.com/vanheeringen-lab/seq2science/compare/v0.3.1...v0.3.0
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
