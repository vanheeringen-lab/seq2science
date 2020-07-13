# Changelog
All notable changes to `seq2science` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

All changed fall under either one of these types: `added`, `changed`, `deprecated`, `removed`, `fixed`, `security`.

## [Unreleased]

### Added

- bwa-mem2 as aligner

### Fixed

- Bug when merging replicates and having controls

## v0.0.3 - 2020-07-01

### Fixed
- bug when specifying 2 cores, which rounded down to zero cores for samtools sorting and crash
- edger environment was incompatible
- seq2science cache on sensible location + seq2science clean fixed

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
