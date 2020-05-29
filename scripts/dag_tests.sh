#!/usr/bin/env bash

# to run these tests locally:
# 1)   chmod u+x ./scripts/dag_tests.sh
# 2)   sh ./scripts/dag_tests.sh

# check if an argument is passed
if [ -z "$1" ]
  then
    echo "No argument supplied"
    exit
fi

CORES=48
trap "rm -rf Jenkins_results; rm -rf Jenkins/dag_fastqs; rm -rf Jenkins/assembly*; rm -rf ~/.config/snakemake/layouts*" EXIT  # remove the test outputs on exit
set -e  # Exit immediately if a command exits with a non-zero status.

if [ $1 = "download_n_alignment" ]; then

  # download workflow
  WF=download_fastq

  printf "\ndownload default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml

  # alignment workflow
  WF=alignment

  mkdir -p Jenkins/dag_fastqs
  touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
  touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz

  printf "\nalignment default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml

  mkdir -p Jenkins/assembly1
  touch Jenkins/assembly1/assembly1.fa
  touch Jenkins/assembly1/assembly1.fa.sizes
  touch Jenkins/assembly1/assembly1.customblacklist_complement.bed
  mkdir -p Jenkins_results/dedup
  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam
  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai

  printf "\ntrackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True

  mkdir -p Jenkins_results/fastq_trimmed
  touch Jenkins_results/fastq_trimmed/S1_1_R1_trimmed.fastq.gz
  touch Jenkins_results/fastq_trimmed/S1_1_R2_trimmed.fastq.gz
  mkdir -p Jenkins_results/qc/trimming
  touch Jenkins_results/qc/trimming/S1_1_R1.fastq.gz_trimming_report.txt
  touch Jenkins_results/qc/trimming/S1_1_R2.fastq.gz_trimming_report.txt

  printf "\nqc multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True

  mkdir -p Jenkins_results/qc/fastqc
  touch Jenkins_results/qc/fastqc/S1_1_R1_fastqc.zip
  touch Jenkins_results/qc/fastqc/S1_1_R2_fastqc.zip
  touch Jenkins_results/qc/fastqc/S1_1_R1_trimmed_fastqc.zip
  touch Jenkins_results/qc/fastqc/S1_1_R2_trimmed_fastqc.zip

  printf "\naligners\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/bowtie2.yaml
  # snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/bwa.yaml  # default
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/hisat2.yaml
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/star.yaml  # will download annotation too

  mkdir -p Jenkins_results/bwa/
  touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam

  printf "\nalignmentsieve\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/alignmentsieve.yaml
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/no_alignmentsieve.yaml

  touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-sieved.bam

  printf "\nsorting\n"
  # snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/samtools_coordinate.yaml  # default
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/samtools_queryname.yaml
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/sambamba_coordinate.yaml
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/sambamba_queryname.yaml

  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam

  printf "\ncram support\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/cram_no_bam.yaml

  touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz  # TODO: bug: these should not be required with the bam file
  touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam

  printf "\nmultiple assemblies\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv

  mkdir -p Jenkins/assembly2
  touch Jenkins/assembly2/assembly2.fa
  touch Jenkins/assembly2/assembly2.fa.sizes
  touch Jenkins/assembly2/assembly2.customblacklist_complement.bed
  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam  # for some reason, this samples's file needs to be recreated(?)
  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai
  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam.bai

  printf "\nmultiple assemblies - trackhubs\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_trackhub=True

  touch Jenkins_results/qc/fastqc/S1_2_R1_fastqc.zip
  touch Jenkins_results/qc/fastqc/S1_2_R2_fastqc.zip
  touch Jenkins_results/qc/fastqc/S1_2_R1_trimmed_fastqc.zip
  touch Jenkins_results/qc/fastqc/S1_2_R2_trimmed_fastqc.zip
  touch Jenkins_results/qc/trimming/S1_2_R1.fastq.gz_trimming_report.txt
  touch Jenkins_results/qc/trimming/S1_2_R2.fastq.gz_trimming_report.txt
  touch Jenkins_results/bwa/assembly2-S1_2.samtools-coordinate-unsieved.bam
  touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
  touch Jenkins_results/bwa/assembly2-S1_2.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam  # for some reason, this samples's file needs to be recreated(?)
  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam.bai
  mkdir -p Jenkins_results/qc/dedup/
  touch Jenkins_results/qc/dedup/assembly1-S1_1.samtools-coordinate.metrics.txt
  touch Jenkins_results/qc/dedup/assembly2-S1_2.samtools-coordinate.metrics.txt

  printf "\nmultiple assemblies - multiqc\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_qc_report=True

  #touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz  # TODO: created above due to bug in layout
  #touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
  touch Jenkins_results/fastq_trimmed/S1_2_R1_trimmed.fastq.gz
  touch Jenkins_results/fastq_trimmed/S1_2_R2_trimmed.fastq.gz

  printf "\nmultiple replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config technical_replicates=merge  # nothing to merge
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=keep
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=merge

  touch Jenkins/dag_fastqs/S2_1.fastq.gz
  touch Jenkins/dag_fastqs/S2_2.fastq.gz

  printf "\nmultiple assemblies and replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=keep

  mkdir -p Jenkins_results/fastq_trimmed/merged
  touch Jenkins_results/fastq_trimmed/S2_1_trimmed.fastq.gz
  touch Jenkins_results/fastq_trimmed/S2_2_trimmed.fastq.gz

  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge

  touch Jenkins_results/dedup/assembly1-R1_1.samtools-coordinate.bam
  touch Jenkins_results/dedup/assembly2-R2_1.samtools-coordinate.bam

  printf "\nmultiple assemblies and replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge create_trackhub=True

  touch Jenkins_results/qc/fastqc/S2_1_fastqc.zip
  touch Jenkins_results/qc/fastqc/S2_2_fastqc.zip
  touch Jenkins_results/qc/trimming/S2_1.fastq.gz_trimming_report.txt
  touch Jenkins_results/qc/trimming/S2_2.fastq.gz_trimming_report.txt
  touch Jenkins_results/bwa/assembly1-R1_1.samtools-coordinate-unsieved.bam
  touch Jenkins_results/bwa/assembly2-R2_1.samtools-coordinate-unsieved.bam
  touch Jenkins_results/dedup/assembly1-R1_1.samtools-coordinate.bam
  touch Jenkins_results/dedup/assembly2-R2_1.samtools-coordinate.bam
  touch Jenkins_results/dedup/assembly1-R1_1.samtools-coordinate.bam.bai
  touch Jenkins_results/dedup/assembly2-R2_1.samtools-coordinate.bam.bai
  touch Jenkins_results/qc/dedup/assembly1-R1_1.samtools-coordinate.metrics.txt
  touch Jenkins_results/qc/dedup/assembly2-R2_1.samtools-coordinate.metrics.txt
  touch Jenkins_results/bwa/assembly1-R1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
  touch Jenkins_results/bwa/assembly2-R2_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge create_qc_report=True

fi

exit
############################################################ <- up untill here seems to work (on my end, with a clean layout)!

## TODO: intermediate files
#touch Jenkins/dag_fastqs/S2_1.fastq.gz
#touch Jenkins/dag_fastqs/S2_2.fastq.gz
#mkdir -p Jenkins_results/fastq_trimmed/merged
#mkdir -p Jenkins_results/qc/samtools_stats
#
#s=2  # sample
#a=2  # assembly
#for r in {1..2}; do  # replicate
##  touch Jenkins_results/qc/fastqc/S${s}_${r}_fastqc.zip
##  touch Jenkins_results/qc/fastqc/S${s}_${r}_trimmed_fastqc.zip
#  touch Jenkins_results/fastq_trimmed/merged/R${s}_${r}_trimmed.fastq.gz
##  touch Jenkins_results/qc/trimming/S${s}_${r}.fastq.gz_trimming_report.txt
#done
#r=1 # replicate
#touch Jenkins_results/dedup/assembly${a}-R${s}_${r}.samtools-coordinate.bam
#touch Jenkins_results/dedup/assembly${a}-R${s}_${r}.samtools-coordinate.bam.bai
#touch Jenkins_results/bwa/assembly${a}-R${s}_${r}.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#touch Jenkins_results/qc/samtools_stats/assembly${a}-R${s}_${r}.samtools-coordinate.samtools_stats.txt
#
#printf "\nmultiple replicates and assemblies\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv create_qc_report=True
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv create_trackhub=True

###################

#mkdir -p Jenkins/assembly2
#touch Jenkins/assembly2/assembly2.fa
#touch Jenkins/assembly2/assembly2.fa.sizes
#touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz
#touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
#touch Jenkins_results/qc/fastqc/S1_2_R1_fastqc.zip
#touch Jenkins_results/qc/fastqc/S1_2_R2_fastqc.zip
#touch Jenkins_results/fastq_trimmed/S1_2_R1_trimmed.fastq.gz
#touch Jenkins_results/fastq_trimmed/S1_2_R2_trimmed.fastq.gz
#touch Jenkins_results/qc/fastqc/S1_2_R1_trimmed_fastqc.zip
#touch Jenkins_results/qc/fastqc/S1_2_R2_trimmed_fastqc.zip
#touch Jenkins_results/qc/trimming/S1_2_R1.fastq.gz_trimming_report.txt
#touch Jenkins_results/qc/trimming/S1_2_R2.fastq.gz_trimming_report.txt

#mkdir -p Jenkins_results/bwa
#touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#mkdir -p Jenkins_results/qc/samtools_stats
#touch Jenkins_results/qc/samtools_stats/assembly1-S1_1.samtools-coordinate.samtools_stats.txt

# --dag | dot -Tsvg > graph.svg
exit

# TODO: should work untill here

#
#
## ATAC-seq workflow
#WF=atac_seq
#
#printf "\ndefault atac-seq\n\n"
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml
#
#printf "\nqc report\n"
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True
#
#printf "\ntrackhub\n"
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True
#
#touch Jenkins_results/dedup/assembly1-S1_1.sambamba-queryname.bam
#touch Jenkins_results/dedup/assembly1-S1_1.sambamba-queryname.bam.bai
#
#printf "\npeak callers\n"
## snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/macs2.yaml  # default
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich.yaml
#
#mkdir -p Jenkins_results/qc/dedup
#touch Jenkins_results/qc/dedup/assembly1-S1_1.samtools-coordinate.metrics.txt
#
#printf "\nmultiple peak callers - qc report\n"
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_qc_report=True
#
#printf "\nmultiple peak callers - trackhub\n"
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_trackhub=True
#
## TODO: intermediate files
### generate data
##mkdir -p Jenkins_results/qc/fastqc
##touch Jenkins_results/qc/fastqc/S1_1_R1_fastqc.zip
##touch Jenkins_results/qc/fastqc/S1_1_R2_fastqc.zip
##touch Jenkins_results/qc/fastqc/S1_1_R1_trimmed_fastqc.zip  # required for minimal run
##touch Jenkins_results/qc/fastqc/S1_1_R2_trimmed_fastqc.zip  # required for minimal run
##mkdir -p Jenkins_results/qc/trimming
##touch Jenkins_results/qc/trimming/S1_1_R1.fastq.gz_trimming_report.txt
##touch Jenkins_results/qc/trimming/S1_1_R2.fastq.gz_trimming_report.txt
##mkdir -p Jenkins_results/bwa
##touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
##mkdir -p Jenkins_results/qc/samtools_stats
##touch Jenkins_results/qc/samtools_stats/assembly1-S1_1.samtools-coordinate.samtools_stats.txt
##mkdir -p Jenkins_results/dedup
##touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam  # required for minimal run
##touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai  # required for minimal run
##mkdir -p Jenkins_results/qc/dedup
##touch Jenkins_results/qc/dedup/assembly1-S1_1.samtools-coordinate.metrics.txt
#
#j=1
#for i in {3..8}; do
#  if (( $i > 4 )); then
#    j=2
#  fi
#  touch Jenkins/dag_fastqs/S${i}_1.fastq.gz
#  touch Jenkins_results/qc/fastqc/S${i}_1_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S${i}_1_trimmed_fastqc.zip
#  touch Jenkins_results/qc/trimming/S${i}_1.fastq.gz_trimming_report.txt
#  touch Jenkins_results/bwa/assembly${j}-S${i}_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#  touch Jenkins_results/dedup/assembly${j}-S${i}_1.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly${j}-S${i}_1.samtools-coordinate.bam.bai
#  touch Jenkins_results/qc/samtools_stats/assembly${j}-S${i}_1.samtools-coordinate.samtools_stats.txt
#  touch Jenkins_results/qc/dedup/assembly${j}-S${i}_1.samtools-coordinate.metrics.txt
#done
#
#printf "\nmultiple replicates and assemblies\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv technical_replicates=keep
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv
#
## TODO: intermediate files
#
#printf "\nmultiple replicates, conditions and assemblies\n"
### snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=keep  # default
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=idr
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=fisher
#
## TODO: intermediate files
#
#printf "\nmultiple replicates, conditions and assemblies - qc report\n"
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=idr create_qc_report=True
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=fisher create_qc_report=True
#
#printf "\nmultiple replicates, conditions and assemblies - trackhub\n"
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=idr create_trackhub=True
##snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv combine_replicates=fisher create_trackhub=True
#
## --dag | dot -Tsvg > graph.svg
