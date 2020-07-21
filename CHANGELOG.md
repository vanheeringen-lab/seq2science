# Changelog
All notable changes to `seq2science` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

All changed fall under either one of these types: `added`, `changed`, `deprecated`, `removed`, `fixed`, `security`.

## [Unreleased]

### Added

- HISAT2 as aligner for RNA-seq
- splice-aware HISAT2 indexing for RNA-seq
- quantifier HTSeq for RNA-seq
- quantifier featurecounts for RNA-seq
- added/expanded `seq2science explain` info (now covers RNA- and scATAC-seq too)

### Changed

- rules and scriptnames in RNA-seq. ex: `txi.R` is now `quant_to_counts.R` to better reflect its function
- `quant_to_counts.R` now converts salmon transcript abundances to gene counts identically to DESeq2
- STAR no longer outputs counts, and is no longer found under `quantifiers`
- gene counts are generated from (filtered) bams when using either STAR or HISAT2 as aligner and HTSeq or featureCounts are quantifier
- batch corrected gene counts are generated if a DESeq2 design contrast inclused a batch
- batch corrected TPM are generated if a DESeq2 design contrast inclused a batch, and quantification was performed using Salmon
  - for us in ANANSE, for instance
- `seq2science explain` now retrieves messages from `explain.smk`.

## v0.1.0 - 2020-07-15

### Added

- bwa-mem2 as aligner
- new command-line option `explain`, which explains what has been done, and writes your material & methods section for you!

### Changed

- change the workflow names, replaced _ by -. (download_fastq to download-fastq, chip_seq to chip-seq, atac_seq to atac-seq, scatac_seq to scatac-seq, and rna_seq to rna-seq)
- changed the way seq2science is called. Moved all the logic from bin/seq2science to seq2science/cli.py

### Fixed

- Bug when merging replicates and having controls

## v0.0.3 - 2020-07-01

### Fixed
- bug when specifying 2 cores, which rounded down to zero cores for samtools sorting and crash
- edger environment was incompatible
- seq2science cache on sensible location + seq2science clean fixed
- only lookup sample layout when not local, opens up for slightly better tests in bioconda recipe

## v0.0.2 - 2020-06-29

### Fixed

- samtools using the correct nr of threads after update to v1.10

### Changed

- The count table for ATAC/ChIP-seq peaks is now made from finding all peaks within a range of 200 bp, and taking the most significant one (gimmemotifs' combine_peaks) and extending the remaining peaks 200 bp. On this count table quantile normalisation, TMM, RLE and upperquartile normalisation with CPM is done. Downstream steps log transform these and mean center them. This however means that for broadpeaks no count_table is generated.
- Snakefmt -l 121 applied

## v0.0.1 - 2020-06-17
Many minor bug- and quality of life fixes.

## v0.0.0 - 2020-06-11
First release of seq2science!
