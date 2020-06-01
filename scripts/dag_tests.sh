#!/usr/bin/env bash

# to run these tests locally:
# 1)   chmod u+x ./scripts/dag_tests.sh
# 2)   ./scripts/dag_tests.sh TEST

if [ -z "$1" ]
  then
    echo "No test specified"
    exit
fi

CORES=48
trap "rm -rf Jenkins_results; rm -rf Jenkins/dag_fastqs; rm -rf ~/.config/snakemake/layouts*" EXIT  # remove the test outputs on exit
set -e  # Exit immediately if a command exits with a non-zero status.
function assert_rulecount {
  # check if the DAG (stored with  | tee Jenkins_results/val  ) ran rule $1 exactly $2 times
  val=$(cat Jenkins_results/val | grep -w $1 | cut -f2);

  # check if the rule was found in the DAG at all
  if [ -z "$val" ]; then
    # if specified count is zero, that's OK
    (($2!=0)) && printf "\nrule $1 did not run at all. Exiting.\n\n" && exit 1;
  else
    # else check if the count equals the specified number
    (($val!=$2)) && printf "\nrule $1 ran $val times instead of the expected $2 times. Exiting.\n\n" && exit 1;
  fi
  :
}
mkdir -p Jenkins_results
mkdir -p Jenkins/dag_fastqs
touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz
touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz
touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
touch Jenkins/dag_fastqs/S2_1.fastq.gz
touch Jenkins/dag_fastqs/S2_2.fastq.gz
touch Jenkins/dag_fastqs/S3_1.fastq.gz
touch Jenkins/dag_fastqs/S4_1.fastq.gz
touch Jenkins/dag_fastqs/S5_1.fastq.gz
touch Jenkins/dag_fastqs/S6_1.fastq.gz
touch Jenkins/dag_fastqs/S7_1.fastq.gz
touch Jenkins/dag_fastqs/S8_1.fastq.gz

if [ $1 = "alignment" ]; then

  # download workflow
  WF=download_fastq

  printf "\ndownload default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml | tee Jenkins_results/val
  assert_rulecount sra2fastq_PE 1
  assert_rulecount sra2fastq_SE 1

  # alignment workflow
  WF=alignment

#  mkdir -p Jenkins/dag_fastqs
#  touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
#  touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz

  printf "\nalignment default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml | tee Jenkins_results/val
  assert_rulecount bwa_index 1
  assert_rulecount samtools_presort 1
  assert_rulecount sambamba_sort 0
  assert_rulecount mark_duplicates 1

#  mkdir -p Jenkins/assembly1
#  touch Jenkins/assembly1/assembly1.fa
#  touch Jenkins/assembly1/assembly1.fa.sizes
#  touch Jenkins/assembly1/assembly1.customblacklist_complement.bed
#  mkdir -p Jenkins_results/dedup
#  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai

  printf "\ntrackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bam_bigwig 1
  assert_rulecount cytoband 1
  assert_rulecount gcPercent 1
  assert_rulecount softmask_track_2 1
  assert_rulecount twobit 1
  assert_rulecount trackhub 1

  printf "\ntrackhub - stranded bams\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/stranded_sample.tsv create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bam_stranded_bigwig 1

#  mkdir -p Jenkins_results/fastq_trimmed
#  touch Jenkins_results/fastq_trimmed/S1_1_R1_trimmed.fastq.gz
#  touch Jenkins_results/fastq_trimmed/S1_1_R2_trimmed.fastq.gz
#  mkdir -p Jenkins_results/qc/trimming
#  touch Jenkins_results/qc/trimming/S1_1_R1.fastq.gz_trimming_report.txt
#  touch Jenkins_results/qc/trimming/S1_1_R2.fastq.gz_trimming_report.txt

  printf "\nmultiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 4
  assert_rulecount multiqc 1

#  mkdir -p Jenkins_results/qc/fastqc
#  touch Jenkins_results/qc/fastqc/S1_1_R1_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_1_R2_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_1_R1_trimmed_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_1_R2_trimmed_fastqc.zip

  printf "\naligners\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/bowtie2.yaml | tee Jenkins_results/val
  assert_rulecount bowtie2_index 1
#  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/bwa.yaml | tee Jenkins_results/val  # default
#  assert_rulecount bwa_index 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/hisat2.yaml | tee Jenkins_results/val
  assert_rulecount hisat2_index 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/star.yaml | tee Jenkins_results/val
  assert_rulecount star_index 1

#  mkdir -p Jenkins_results/bwa/
#  touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam

  printf "\nalignmentsieve\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/alignmentsieve.yaml | tee Jenkins_results/val
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/no_alignmentsieve.yaml | tee Jenkins_results/val

#  touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-sieved.bam

  printf "\nsorting\n"
  # snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/samtools_coordinate.yaml  # default
  # assert_rulecount samtools_presort 1
  # assert_rulecount sambamba_sort 0
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/samtools_queryname.yaml | tee Jenkins_results/val
  assert_rulecount samtools_sort 1
  assert_rulecount sambamba_sort 0
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/sambamba_coordinate.yaml | tee Jenkins_results/val
  assert_rulecount sambamba_sort 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/sambamba_queryname.yaml | tee Jenkins_results/val
  assert_rulecount sambamba_sort 1

#  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam

  printf "\ncram support\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/cram_no_bam.yaml | tee Jenkins_results/val
  assert_rulecount bam2cram 1

#  touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz
#  touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam

  printf "\nmultiple assemblies\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv | tee Jenkins_results/val
  assert_rulecount bwa_index 2
  assert_rulecount bwa_mem 2

#  mkdir -p Jenkins/assembly2
#  touch Jenkins/assembly2/assembly2.fa
#  touch Jenkins/assembly2/assembly2.fa.sizes
#  touch Jenkins/assembly2/assembly2.customblacklist_complement.bed
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam  # for some reason, this samples's file needs to be recreated(?)
#  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam.bai
#  touch Jenkins/assembly1/cytoBandIdeo.bb
#  touch Jenkins/assembly2/cytoBandIdeo.bb
#  touch Jenkins/assembly1/assembly1.gc5Base.bw
#  touch Jenkins/assembly2/assembly2.gc5Base.bw
#  touch Jenkins/assembly1/assembly1_softmasking.bb
#  touch Jenkins/assembly2/assembly2_softmasking.bb

  printf "\nmultiple assemblies - trackhubs\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_trackhub=True | tee Jenkins_results/val
  assert_rulecount twobit 2
  assert_rulecount trackhub 1

#  touch Jenkins_results/qc/fastqc/S1_2_R1_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_2_R2_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_2_R1_trimmed_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_2_R2_trimmed_fastqc.zip
#  touch Jenkins_results/qc/trimming/S1_2_R1.fastq.gz_trimming_report.txt
#  touch Jenkins_results/qc/trimming/S1_2_R2.fastq.gz_trimming_report.txt
#  touch Jenkins_results/bwa/assembly2-S1_2.samtools-coordinate-unsieved.bam
#  touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#  touch Jenkins_results/bwa/assembly2-S1_2.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam  # for some reason, this samples's file needs to be recreated(?)
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam.bai
#  mkdir -p Jenkins_results/qc/dedup/
#  touch Jenkins_results/qc/dedup/assembly1-S1_1.samtools-coordinate.metrics.txt
#  touch Jenkins_results/qc/dedup/assembly2-S1_2.samtools-coordinate.metrics.txt

  printf "\nmultiple assemblies - multiqc\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_qc_report=True | tee Jenkins_results/val
  assert_rulecount multiqc 2

#  touch Jenkins_results/fastq_trimmed/S1_2_R1_trimmed.fastq.gz
#  touch Jenkins_results/fastq_trimmed/S1_2_R2_trimmed.fastq.gz

  printf "\nmultiple replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config technical_replicates=merge | tee Jenkins_results/val  # nothing to merge
  assert_rulecount merge_replicates 0
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount merge_replicates 0
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=merge | tee Jenkins_results/val
  assert_rulecount merge_replicates 2
  assert_rulecount bwa_mem 1

#  touch Jenkins/dag_fastqs/S2_1.fastq.gz
#  touch Jenkins/dag_fastqs/S2_2.fastq.gz

  printf "\nmultiple assemblies and replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount bwa_index 2
  assert_rulecount bwa_mem 4

#  mkdir -p Jenkins_results/fastq_trimmed/merged
#  touch Jenkins_results/fastq_trimmed/S2_1_trimmed.fastq.gz
#  touch Jenkins_results/fastq_trimmed/S2_2_trimmed.fastq.gz

  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge | tee Jenkins_results/val
  assert_rulecount bwa_index 2
  assert_rulecount bwa_mem 2

#  touch Jenkins_results/dedup/assembly1-R1_1.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly2-R2_1.samtools-coordinate.bam

  printf "\nmultiple assemblies and replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bwa_mem 2
  assert_rulecount bam_bigwig 2
#
#  touch Jenkins_results/qc/fastqc/S2_1_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S2_2_fastqc.zip
#  touch Jenkins_results/qc/trimming/S2_1.fastq.gz_trimming_report.txt
#  touch Jenkins_results/qc/trimming/S2_2.fastq.gz_trimming_report.txt
#  touch Jenkins_results/bwa/assembly1-R1_1.samtools-coordinate-unsieved.bam
#  touch Jenkins_results/bwa/assembly2-R2_1.samtools-coordinate-unsieved.bam
#  touch Jenkins_results/dedup/assembly1-R1_1.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly2-R2_1.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly1-R1_1.samtools-coordinate.bam.bai
#  touch Jenkins_results/dedup/assembly2-R2_1.samtools-coordinate.bam.bai
#  touch Jenkins_results/qc/dedup/assembly1-R1_1.samtools-coordinate.metrics.txt
#  touch Jenkins_results/qc/dedup/assembly2-R2_1.samtools-coordinate.metrics.txt
#  touch Jenkins_results/bwa/assembly1-R1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#  touch Jenkins_results/bwa/assembly2-R2_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
  assert_rulecount samtools_stats 2

fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq

#  mkdir -p Jenkins/dag_fastqs
#  touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz  # layout dependency
#  touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz
#
#  mkdir -p Jenkins/assembly1
#  touch Jenkins/assembly1/assembly1.fa
#  touch Jenkins/assembly1/assembly1.fa.sizes
#  touch Jenkins/assembly1/assembly1.customblacklist_complement.bed
#  touch Jenkins/assembly1/assembly1.2bit
#  touch Jenkins/assembly1/cytoBandIdeo.bb
#  touch Jenkins/assembly1/assembly1.gc5Base.bw
#  touch Jenkins/assembly1/assembly1_softmasking.bb
#
#  mkdir -p Jenkins_results/dedup
#  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai
#
#  mkdir -p Jenkins_results/qc/fastqc
#  touch Jenkins_results/qc/fastqc/S1_1_R1_trimmed_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_1_R2_trimmed_fastqc.zip

  printf "\natac-seq default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml

#  mkdir -p Jenkins_results/macs2
#  touch Jenkins_results/macs2/assembly1-S1_1_control_lambda.bdg
#  touch Jenkins_results/macs2/assembly1-S1_1_peaks.xls
#  touch Jenkins_results/macs2/assembly1-S1_1_treat_pileup.bdg
#  touch Jenkins_results/macs2/assembly1-S1_1_summits.bed
#  touch Jenkins_results/macs2/assembly1-S1_1_peaks.narrowPeak
#  touch Jenkins_results/macs2/assembly1_peaks.bed
#
#  mkdir -p Jenkins_results/count_table/macs2/
#  touch Jenkins_results/count_table/macs2/count_table_assembly1.samtools-coordinate.txt

  printf "\ntrackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True

#  touch Jenkins_results/qc/fastqc/S1_1_R1_fastqc.zip
#  touch Jenkins_results/qc/fastqc/S1_1_R2_fastqc.zip
#  mkdir -p Jenkins_results/qc/trimming
#  touch Jenkins_results/qc/trimming/S1_1_R1.fastq.gz_trimming_report.txt
#  touch Jenkins_results/qc/trimming/S1_1_R2.fastq.gz_trimming_report.txt
#
#  mkdir -p touch Jenkins_results/qc/samtools_stats
#  touch Jenkins_results/qc/samtools_stats/assembly1-S1_1.samtools-coordinate.samtools_stats.txt
#
#  mkdir -p Jenkins_results/bwa
#  touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#
#  mkdir -p Jenkins_results/qc/dedup
#  touch Jenkins_results/qc/dedup/assembly1-S1_1.samtools-coordinate.metrics.txt

  printf "\nqc multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True

#  touch Jenkins_results/dedup/assembly1-S1_1.sambamba-queryname.bam
#  touch Jenkins_results/dedup/assembly1-S1_1.sambamba-queryname.bam.bai

  printf "\npeak callers\n"
  # snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/macs2.yaml  # default
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich.yaml

#  rm Jenkins_results/count_table/macs2/count_table_assembly1.samtools-coordinate.txt

  printf "\nmultiple peak callers\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml

#  mkdir -p Jenkins_results/genrich
#  touch Jenkins_results/genrich/assembly1-S1_1.log
#  touch Jenkins_results/genrich/assembly1-S1_1_peaks.narrowPeak
#
#  mkdir -p Jenkins_results/count_table/genrich
#  touch Jenkins_results/count_table/macs2/count_table_assembly1.samtools-coordinate.txt
#  touch Jenkins_results/count_table/genrich/count_table_assembly1.samtools-coordinate.txt

  printf "\nmultiple peak callers - qc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_qc_report=True

#  touch Jenkins_results/genrich/assembly1-S1_1.bedgraph

  printf "\nmultiple peak callers - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_trackhub=True

#  touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz
#  touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam.bai
#  touch Jenkins_results/macs2/assembly2_peaks.bed
#  touch Jenkins_results/genrich/assembly2_peaks.bed

  printf "\nmultiple peak callers & multiple assemblies\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/alignment/assemblies.tsv

#  mkdir -p Jenkins/assembly2
#  touch Jenkins/assembly2/assembly2.fa
#  touch Jenkins/assembly2/assembly2.fa.sizes
#  touch Jenkins/assembly2/assembly2.customblacklist_complement.bed
#  touch Jenkins/assembly2/assembly2.2bit
#  touch Jenkins/assembly2/cytoBandIdeo.bb
#  touch Jenkins/assembly2/assembly2.gc5Base.bw
#  touch Jenkins/assembly2/assembly2_softmasking.bb
#
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam
#  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam.bai
#
#  touch Jenkins_results/genrich/assembly2-S1_2.log
#  touch Jenkins_results/genrich/assembly2-S1_2.bdgish
#  touch Jenkins_results/genrich/assembly2-S1_2_peaks.narrowPeak
#  touch Jenkins_results/macs2/assembly2-S1_2_control_lambda.bdg
#  touch Jenkins_results/macs2/assembly2-S1_2_peaks.narrowPeak
#  touch Jenkins_results/macs2/assembly2-S1_2_treat_pileup.bdg
#
#  touch Jenkins_results/macs2/assembly2_peaks.bed
#  touch Jenkins_results/genrich/assembly2_peaks.bed
#
#  touch Jenkins_results/count_table/genrich/count_table_assembly2.samtools-coordinate.txt
#  touch Jenkins_results/count_table/macs2/count_table_assembly2.samtools-coordinate.txt

  printf "\nmultiple peak callers & multiple assemblies - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_trackhub=True

# incomplete intermediate files
#touch Jenkins_results/qc/fastqc/S1_2_fastqc.zip
#touch Jenkins_results/qc/fastqc/S1_2_fastqc.zip
#touch Jenkins_results/qc/trimming/S1_2.fastq.gz_trimming_report.txt
#touch Jenkins_results/qc/trimming/S1_2.fastq.gz_trimming_report.txt
##touch Jenkins_results/qc/samtools_stats/assembly2-S1_2.samtools-coordinate.samtools_stats.txt
#touch Jenkins_results/bwa/assembly2-S1_2.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json
#touch Jenkins_results/qc/dedup/assembly2-S1_2.samtools-coordinate.metrics.txt

  printf "\nmultiple peak callers & multiple assemblies - qc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_qc_report=True

  printf "\nmultiple peak callers & multiple replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/alignment/replicates.tsv

  printf "\nmultiple peak callers & multiple replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/alignment/replicates.tsv create_trackhub=True

  printf "\nmultiple peak callers & multiple replicates - qc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/alignment/replicates.tsv create_qc_report=True

#  touch Jenkins/dag_fastqs/S2_1.fastq.gz
#  touch Jenkins/dag_fastqs/S2_2.fastq.gz
#  touch Jenkins/dag_fastqs/S3_1.fastq.gz
#  touch Jenkins/dag_fastqs/S4_1.fastq.gz
#  touch Jenkins/dag_fastqs/S5_1.fastq.gz
#  touch Jenkins/dag_fastqs/S6_1.fastq.gz
#  touch Jenkins/dag_fastqs/S7_1.fastq.gz
#  touch Jenkins/dag_fastqs/S8_1.fastq.gz

  printf "\nmultiple peak callers, assemblies and replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv

  printf "\nmultiple peak callers, assemblies and replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv create_trackhub=True

  printf "\nmultiple peak callers, assemblies and replicates - qc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv create_qc_report=True

fi

if [ $1 = "scatac-seq" ]; then

# ATAC-seq workflow (also covers ChIP-seq workflow)
WF=scATAC_seq

#mkdir -p Jenkins_results
#mkdir -p Jenkins/dag_fastqs
#touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
#touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz
#touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz
#touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
#touch Jenkins/dag_fastqs/S2_1.fastq.gz
#touch Jenkins/dag_fastqs/S2_2.fastq.gz
#touch Jenkins/dag_fastqs/S3_1.fastq.gz
#touch Jenkins/dag_fastqs/S4_1.fastq.gz
#touch Jenkins/dag_fastqs/S5_1.fastq.gz
#touch Jenkins/dag_fastqs/S6_1.fastq.gz
#touch Jenkins/dag_fastqs/S7_1.fastq.gz
#touch Jenkins/dag_fastqs/S8_1.fastq.gz

printf "\nscatac-seq default\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml | tee Jenkins_results/val
assert_rulecount create_SNAP_object 1

printf "\ntrackhub\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config create_trackhub=True | tee Jenkins_results/val
# TODO: scATAC-seq does not create a trackhub
#assert_rulecount trackhub 1

printf "\nqc multiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config create_qc_report=True

printf "\nmultiple assemblies\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv

printf "\nmultiple assemblies - trackhubs\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_trackhub=True

printf "\nmultiple assemblies - multiqc\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_qc_report=True

printf "\nmultiple replicates\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config technical_replicates=merge  # nothing to merge
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=keep
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=merge

printf "\nmultiple replicates - trackhub\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=merge create_trackhub=True

printf "\nmultiple replicates - multiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=merge create_qc_report=True

printf "\nmultiple assemblies and replicates\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=keep
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge

printf "\nmultiple assemblies and replicates - trackhub\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge create_trackhub=True

printf "\nmultiple assemblies and replicates - multiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge create_qc_report=True

fi

if [ $1 = "rna-seq" ]; then

# RNA-seq workflow
WF=rna_seq

#mkdir -p Jenkins/dag_fastqs
#touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
#touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz

printf "\nrna-seq default\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=star

printf "\nquantifiers\n"
# snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=star  # default
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=salmon

printf "\ntrackhub\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=star create_trackhub=True
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=salmon create_trackhub=True

printf "\nqc multiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=star create_qc_report=True

printf "\ndecoy aware salmon index\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=salmon decoy_aware_index=True

#touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz
#touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
#touch Jenkins/dag_fastqs/S2_1.fastq.gz
#touch Jenkins/dag_fastqs/S2_2.fastq.gz
#touch Jenkins/dag_fastqs/S3_1.fastq.gz
#touch Jenkins/dag_fastqs/S4_1.fastq.gz
#touch Jenkins/dag_fastqs/S5_1.fastq.gz
#touch Jenkins/dag_fastqs/S6_1.fastq.gz
#touch Jenkins/dag_fastqs/S7_1.fastq.gz
#touch Jenkins/dag_fastqs/S8_1.fastq.gz

printf "\ndifferential expression analysis\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config quantifier=star technical_replicates=keep
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config quantifier=salmon technical_replicates=keep

printf "\nmultiple assemblies with DEA\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep

printf "\nmultiple assemblies with DEA - trackhubs\n"
# TODO: error!
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep create_trackhub=True

printf "\nmultiple assemblies with DEA - multiqc\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep create_qc_report=True

printf "\nmultiple replicates with DEA \n"
# snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config technical_replicates=keep  # default
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config technical_replicates=merge

printf "\nmultiple assemblies and replicates with DEA \n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=keep quantifier=star
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=star
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=salmon

printf "\nmultiple assemblies and replicates with DEA - trackhub\n"
# TODO: error!
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=merge create_trackhub=True

printf "\nmultiple assemblies and replicates with DEA - multiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=merge create_qc_report=True

fi
