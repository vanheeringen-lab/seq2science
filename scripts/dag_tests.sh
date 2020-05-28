#!/usr/bin/env bash

# to run these tests locally:
# 1)   chmod u+x ./scripts/dag_tests.sh
# 2)   ./scripts/dag_tests.sh

CORES=48
trap "rm -rf Jenkins_results; rm -rf Jenkins/dag_fastqs; rm -rf Jenkins/assembly1; rm -rf Jenkins/assembly2" EXIT  # remove the test outputs on exit
set -e  # Exit immediately if a command exits with a non-zero status.

## download workflow
#WF=download_fastq
#printf "\ndownload default\n\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml

# alignment workflow
WF=alignment

# TODO: test genome & annotation downloading in DAG > requires extra samples.tsv
mkdir -p Jenkins/assembly1
touch Jenkins/assembly1/assembly1.fa
touch Jenkins/assembly1/assembly1.fa.sizes
mkdir -p Jenkins/assembly2
touch Jenkins/assembly2/assembly2.fa
touch Jenkins/assembly2/assembly2.fa.sizes
mkdir -p Jenkins/dag_fastqs
touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz

printf "\nalignment default\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml

mkdir -p Jenkins_results/fastq_trimmed
touch Jenkins_results/fastq_trimmed/S1_1_R1_trimmed.fastq.gz
touch Jenkins_results/fastq_trimmed/S1_1_R2_trimmed.fastq.gz
mkdir -p Jenkins_results/qc/trimming
touch Jenkins_results/qc/trimming/S1_1_R1.fastq.gz_trimming_report.txt
touch Jenkins_results/qc/trimming/S1_1_R2.fastq.gz_trimming_report.txt

printf "\nQC report\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True

mkdir -p Jenkins_results/qc/fastqc
touch Jenkins_results/qc/fastqc/S1_1_R1_fastqc.zip
touch Jenkins_results/qc/fastqc/S1_1_R2_fastqc.zip
touch Jenkins_results/qc/fastqc/S1_1_R1_trimmed_fastqc.zip
touch Jenkins_results/qc/fastqc/S1_1_R2_trimmed_fastqc.zip
mkdir -p Jenkins_results/dedup
touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam
touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai
mkdir -p Jenkins_results/bwa
touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
mkdir -p Jenkins_results/qc/samtools_stats
touch Jenkins_results/qc/samtools_stats/assembly1-S1_1.samtools-coordinate.samtools_stats.txt

printf "\ntrackhub\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True

printf "\naligners\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/bowtie2.yaml
## snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/bwa.yaml  # default
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/hisat2.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/star.yaml

printf "\nsorting\n"
## snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/samtools_coordinate.yaml  # default
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/samtools_queryname.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/sambamba_coordinate.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/sambamba_queryname.yaml

printf "\nalignmentsieve\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/alignmentsieve.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/no_alignmentsieve.yaml

printf "\ncram support\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/cram_no_bam.yaml

touch Jenkins/dag_fastqs/S1_2.fastq.gz      # can we merge a PE and an SE sample?
#touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz  # can we merge a PE and an SE sample?
#touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz  # can we merge a PE and an SE sample?
touch Jenkins/dag_fastqs/S2_1.fastq.gz
touch Jenkins/dag_fastqs/S2_2.fastq.gz

printf "\nreplicates\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=keep
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge

s=2  # sample
a=2  # assembly
for r in {1..2}; do  # replicate
  touch Jenkins_results/qc/fastqc/S${s}_${r}_fastqc.zip
  touch Jenkins_results/qc/fastqc/S${s}_${r}_trimmed_fastqc.zip
  touch Jenkins_results/fastq_trimmed/S${s}_${r}_trimmed.fastq.gz
  touch Jenkins_results/qc/trimming/S${s}_${r}.fastq.gz_trimming_report.txt
done
r=1 # replicate
touch Jenkins_results/dedup/assembly${a}-S${s}_${r}.samtools-coordinate.bam
touch Jenkins_results/dedup/assembly${a}-S${s}_${r}.samtools-coordinate.bam.bai
touch Jenkins_results/bwa/assembly${a}-S${s}_${r}.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
touch Jenkins_results/qc/samtools_stats/assembly${a}-S${s}_${r}.samtools-coordinate.samtools_stats.txt

printf "\nmultiple assemblies\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv create_qc_report=True
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv create_trackhub=True



#  touch Jenkins_results/fastq_trimmed/S1_${r}_R1_trimmed.fastq.gz
#  touch Jenkins_results/fastq_trimmed/S1_${r}_R2_trimmed.fastq.gz
#  touch Jenkins_results/qc/trimming/S1_${r}_R1.fastq.gz_trimming_report.txt
#  touch Jenkins_results/qc/trimming/S1_${r}_R2.fastq.gz_trimming_report.txt
#  touch Jenkins_results/qc/fastqc/S1_${r}_R1_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_${r}_R2_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_${r}_R1_trimmed_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_${r}_R2_trimmed_fastqc.zip
#  touch Jenkins_results/dedup/assembly${a}-S1_${r}.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly${a}-S1_${r}.samtools-coordinate.bam.bai
#  touch Jenkins_results/bwa/assembly${a}-S1_${r}.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#  touch Jenkins_results/qc/samtools_stats/assembly${a}-S1_${r}.samtools-coordinate.samtools_stats.txt

#i=1
#j=1
#for k in {1..2}; do
#  touch Jenkins/dag_fastqs/S${i}_1_R${$k}.fastq.gz
#  touch Jenkins_results/qc/fastqc/S${i}_1_R${$k}_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S${i}_1_R${$k}_trimmed_fastqc.zip
#  touch Jenkins_results/qc/trimming/S${i}_1_R${$k}.fastq.gz_trimming_report.txt
#  touch Jenkins_results/bwa/assembly${j}-S${i}_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#  touch Jenkins_results/dedup/assembly${j}-S${i}_1.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly${j}-S${i}_1.samtools-coordinate.bam.bai
#  touch Jenkins_results/qc/samtools_stats/assembly${j}-S${i}_1.samtools-coordinate.samtools_stats.txt
#  touch Jenkins_results/qc/dedup/assembly${j}-S${i}_1.samtools-coordinate.metrics.txt
#done

# generate data
mkdir -p Jenkins_results/qc/fastqc
touch Jenkins_results/qc/fastqc/S1_1_R1_fastqc.zip
touch Jenkins_results/qc/fastqc/S1_1_R2_fastqc.zip
touch Jenkins_results/qc/fastqc/S1_1_R1_trimmed_fastqc.zip  # required for minimal run
touch Jenkins_results/qc/fastqc/S1_1_R2_trimmed_fastqc.zip  # required for minimal run
mkdir -p Jenkins_results/qc/trimming
touch Jenkins_results/qc/trimming/S1_1_R1.fastq.gz_trimming_report.txt
touch Jenkins_results/qc/trimming/S1_1_R2.fastq.gz_trimming_report.txt
mkdir -p Jenkins_results/bwa
touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
mkdir -p Jenkins_results/qc/samtools_stats
touch Jenkins_results/qc/samtools_stats/assembly1-S1_1.samtools-coordinate.samtools_stats.txt
mkdir -p Jenkins_results/dedup
touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam  # required for minimal run
touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai  # required for minimal run
mkdir -p Jenkins_results/qc/dedup
touch Jenkins_results/qc/dedup/assembly1-S1_1.samtools-coordinate.metrics.txt

# ATAC-seq workflow
WF=atac_seq
#printf "\ndefault atac-seq\n\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True
#
#printf "\ntest peakcallers\n\n"
#touch Jenkins_results/dedup/assembly1-S1_1.sambamba-queryname.bam
#touch Jenkins_results/dedup/assembly1-S1_1.sambamba-queryname.bam.bai
## snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/macs2.yaml  # default
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_qc_report=True
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_trackhub=True

j=1
for i in {3..8}; do
  if (( $i > 4 )); then
    j=2
  fi
  touch Jenkins/dag_fastqs/S${i}_1.fastq.gz
  touch Jenkins_results/qc/fastqc/S${i}_1_fastqc.zip
  touch Jenkins_results/qc/fastqc/S${i}_1_trimmed_fastqc.zip
  touch Jenkins_results/qc/trimming/S${i}_1.fastq.gz_trimming_report.txt
  touch Jenkins_results/bwa/assembly${j}-S${i}_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
  touch Jenkins_results/dedup/assembly${j}-S${i}_1.samtools-coordinate.bam
  touch Jenkins_results/dedup/assembly${j}-S${i}_1.samtools-coordinate.bam.bai
  touch Jenkins_results/qc/samtools_stats/assembly${j}-S${i}_1.samtools-coordinate.samtools_stats.txt
  touch Jenkins_results/qc/dedup/assembly${j}-S${i}_1.samtools-coordinate.metrics.txt
done

printf "\natac-seq complex\n\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv technical_replicates=keep
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv
## snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=keep  # default
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=idr
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=fisher
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv create_qc_report=True
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv create_trackhub=True

# --dag | dot -Tsvg > graph.svg
