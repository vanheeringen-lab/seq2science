#!/usr/bin/env bash

# to run these tests locally:
# 1)   chmod u+x ./scripts/dag_tests.sh
# 2)   ./scripts/dag_tests.sh

CORES=48
trap "rm -rf Jenkins_results; rm -rf Jenkins/tinydata2" EXIT  # remove the test outputs on exit
set -e  # Exit immediately if a command exits with a non-zero status.

## download workflow
#WF=download_fastq
#printf "\ndownload default\n\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml
#
#
## alignment workflow
#WF=alignment
#printf "\nalignment default\n\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True
#
#printf "\ntest aligners\n\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/bowtie2.yaml
## snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/bwa.yaml  # default
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/hisat2.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/star.yaml
#
#printf "\ntest sorting\n\n"
## snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/samtools_coordinate.yaml  # default
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/samtools_queryname.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/sambamba_coordinate.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/sambamba_queryname.yaml
#
#printf "\ntest alignmentsieve\n\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/alignmentsieve.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/no_alignmentsieve.yaml
#
#printf "\ntest cram support\n\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/cram_no_bam.yaml
#
#printf "\nalignment complex\n\n"
#mkdir Jenkins/tinydata2
#touch Jenkins/tinydata2/tinydata2.fa
#touch Jenkins/tinydata2/tinydata2.fa.sizes
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=keep
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv create_qc_report=True
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv create_trackhub=True


# ATAC-seq workflow
WF=atac_seq
printf "\ndefault atac-seq\n\n"
mkdir -p Jenkins_results/qc/fastqc
touch Jenkins_results/qc/fastqc/T1_1_R1_fastqc.zip
touch Jenkins_results/qc/fastqc/T1_1_R2_fastqc.zip
touch Jenkins_results/qc/fastqc/T1_1_R1_trimmed_fastqc.zip  # required for minimal run
touch Jenkins_results/qc/fastqc/T1_1_R2_trimmed_fastqc.zip  # required for minimal run
mkdir -p Jenkins_results/qc/trimming
touch Jenkins_results/qc/trimming/T1_1_R1.fastq.gz_trimming_report.txt
touch Jenkins_results/qc/trimming/T1_1_R2.fastq.gz_trimming_report.txt
mkdir -p Jenkins_results/bwa
touch Jenkins_results/bwa/tinydata-T1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
mkdir -p Jenkins_results/qc/samtools_stats
touch Jenkins_results/qc/samtools_stats/tinydata-T1_1.samtools-coordinate.samtools_stats.txt
mkdir -p Jenkins_results/dedup
touch Jenkins_results/dedup/tinydata-T1_1.samtools-coordinate.bam  # required for minimal run
touch Jenkins_results/dedup/tinydata-T1_1.samtools-coordinate.bam.bai  # required for minimal run
mkdir -p Jenkins_results/qc/dedup
touch Jenkins_results/qc/dedup/tinydata-T1_1.samtools-coordinate.metrics.txt
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True

printf "\ntest peakcallers\n\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/macs2.yaml  # default
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich.yaml
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_qc_report=True
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_trackhub=True

#printf "\nalignment complex\n\n"
# replicates, conditions, descriptive names
# replicates keep/merge
# condition combining
# with partial descriptive names



## test different replicate settings
#printf "\natac combine replicates\n\n"
#snakemake -s workflows/atac_seq/Snakefile --directory workflows/atac_seq -n -j 48 --config combine_replicates=fisher samples=../../Jenkins/atac_seq/samples.tsv --quiet
#snakemake -s workflows/atac_seq/Snakefile --directory workflows/atac_seq -n -j 48 --config combine_replicates=idr    samples=../../Jenkins/atac_seq/samples.tsv --quiet
#snakemake -s workflows/atac_seq/Snakefile --directory workflows/atac_seq -n -j 48 --config combine_replicates=merge  samples=../../Jenkins/atac_seq/samples.tsv --quiet
#