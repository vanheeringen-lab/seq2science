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
trap "rm -rf Jenkins_results; rm -rf Jenkins/dag_fastqs; rm -rf Jenkins/assembly*; rm -rf ~/.config/snakemake/layouts*" EXIT  # remove the test outputs on exit
set -e  # Exit immediately if a command exits with a non-zero status.

if [ $1 = "alignment" ]; then

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

  printf "\ntrackhub - stranded bams\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config samples=../../Jenkins/alignment/stranded_sample.tsv create_trackhub=True

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
  touch Jenkins/assembly1/cytoBandIdeo.bb
  touch Jenkins/assembly2/cytoBandIdeo.bb
  touch Jenkins/assembly1/assembly1.gc5Base.bw
  touch Jenkins/assembly2/assembly2.gc5Base.bw
  touch Jenkins/assembly1/assembly1_softmasking.bb
  touch Jenkins/assembly2/assembly2_softmasking.bb

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

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq

  mkdir -p Jenkins/dag_fastqs
  touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz  # layout dependency
  touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz

  mkdir -p Jenkins/assembly1
  touch Jenkins/assembly1/assembly1.fa
  touch Jenkins/assembly1/assembly1.fa.sizes
  touch Jenkins/assembly1/assembly1.customblacklist_complement.bed
  touch Jenkins/assembly1/assembly1.2bit
  touch Jenkins/assembly1/cytoBandIdeo.bb
  touch Jenkins/assembly1/assembly1.gc5Base.bw
  touch Jenkins/assembly1/assembly1_softmasking.bb

  mkdir -p Jenkins_results/dedup
  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam
  touch Jenkins_results/dedup/assembly1-S1_1.samtools-coordinate.bam.bai

  mkdir -p Jenkins_results/qc/fastqc
  touch Jenkins_results/qc/fastqc/S1_1_R1_trimmed_fastqc.zip
  touch Jenkins_results/qc/fastqc/S1_1_R2_trimmed_fastqc.zip

  printf "\natac-seq default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml

  mkdir -p Jenkins_results/macs2
  touch Jenkins_results/macs2/assembly1-S1_1_control_lambda.bdg
  touch Jenkins_results/macs2/assembly1-S1_1_peaks.xls
  touch Jenkins_results/macs2/assembly1-S1_1_treat_pileup.bdg
  touch Jenkins_results/macs2/assembly1-S1_1_summits.bed
  touch Jenkins_results/macs2/assembly1-S1_1_peaks.narrowPeak
  touch Jenkins_results/macs2/assembly1_peaks.bed

  mkdir -p Jenkins_results/count_table/macs2/
  touch Jenkins_results/count_table/macs2/count_table_assembly1.samtools-coordinate.txt

  printf "\ntrackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_trackhub=True

  touch Jenkins_results/qc/fastqc/S1_1_R1_fastqc.zip
  touch Jenkins_results/qc/fastqc/S1_1_R2_fastqc.zip
  mkdir -p Jenkins_results/qc/trimming
  touch Jenkins_results/qc/trimming/S1_1_R1.fastq.gz_trimming_report.txt
  touch Jenkins_results/qc/trimming/S1_1_R2.fastq.gz_trimming_report.txt

  mkdir -p touch Jenkins_results/qc/samtools_stats
  touch Jenkins_results/qc/samtools_stats/assembly1-S1_1.samtools-coordinate.samtools_stats.txt

  mkdir -p Jenkins_results/bwa
  touch Jenkins_results/bwa/assembly1-S1_1.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json

  mkdir -p Jenkins_results/qc/dedup
  touch Jenkins_results/qc/dedup/assembly1-S1_1.samtools-coordinate.metrics.txt

  printf "\nqc multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/default_config.yaml --config create_qc_report=True

  touch Jenkins_results/dedup/assembly1-S1_1.sambamba-queryname.bam
  touch Jenkins_results/dedup/assembly1-S1_1.sambamba-queryname.bam.bai

  printf "\npeak callers\n"
  # snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/macs2.yaml  # default
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich.yaml

  rm Jenkins_results/count_table/macs2/count_table_assembly1.samtools-coordinate.txt

  printf "\nmultiple peak callers\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml

  mkdir -p Jenkins_results/genrich
  touch Jenkins_results/genrich/assembly1-S1_1.log
  touch Jenkins_results/genrich/assembly1-S1_1_peaks.narrowPeak

  mkdir -p Jenkins_results/count_table/genrich
  touch Jenkins_results/count_table/macs2/count_table_assembly1.samtools-coordinate.txt
  touch Jenkins_results/count_table/genrich/count_table_assembly1.samtools-coordinate.txt

  printf "\nmultiple peak callers - qc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_qc_report=True

  touch Jenkins_results/genrich/assembly1-S1_1.bedgraph

  printf "\nmultiple peak callers - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config create_trackhub=True

  touch Jenkins/dag_fastqs/S1_2_R1.fastq.gz
  touch Jenkins/dag_fastqs/S1_2_R2.fastq.gz
  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam
  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam.bai
  touch Jenkins_results/macs2/assembly2_peaks.bed
  touch Jenkins_results/genrich/assembly2_peaks.bed

  printf "\nmultiple peak callers & multiple assemblies\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/alignment/assemblies.tsv

  mkdir -p Jenkins/assembly2
  touch Jenkins/assembly2/assembly2.fa
  touch Jenkins/assembly2/assembly2.fa.sizes
  touch Jenkins/assembly2/assembly2.customblacklist_complement.bed
  touch Jenkins/assembly2/assembly2.2bit
  touch Jenkins/assembly2/cytoBandIdeo.bb
  touch Jenkins/assembly2/assembly2.gc5Base.bw
  touch Jenkins/assembly2/assembly2_softmasking.bb

  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam
  touch Jenkins_results/dedup/assembly2-S1_2.samtools-coordinate.bam.bai

  touch Jenkins_results/genrich/assembly2-S1_2.log
  touch Jenkins_results/genrich/assembly2-S1_2.bdgish
  touch Jenkins_results/genrich/assembly2-S1_2_peaks.narrowPeak
  touch Jenkins_results/macs2/assembly2-S1_2_control_lambda.bdg
  touch Jenkins_results/macs2/assembly2-S1_2_peaks.narrowPeak
  touch Jenkins_results/macs2/assembly2-S1_2_treat_pileup.bdg

  touch Jenkins_results/macs2/assembly2_peaks.bed
  touch Jenkins_results/genrich/assembly2_peaks.bed

  touch Jenkins_results/count_table/genrich/count_table_assembly2.samtools-coordinate.txt
  touch Jenkins_results/count_table/macs2/count_table_assembly2.samtools-coordinate.txt

  printf "\nmultiple peak callers & multiple assemblies - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_trackhub=True

  # TODO; no intermediate files added from here. would create clearer DAGs and maybe speed up the test, but its a lot of work :-/
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

  touch Jenkins/dag_fastqs/S2_1.fastq.gz
  touch Jenkins/dag_fastqs/S2_2.fastq.gz
  touch Jenkins/dag_fastqs/S3_1.fastq.gz
  touch Jenkins/dag_fastqs/S4_1.fastq.gz
  touch Jenkins/dag_fastqs/S5_1.fastq.gz
  touch Jenkins/dag_fastqs/S6_1.fastq.gz
  touch Jenkins/dag_fastqs/S7_1.fastq.gz
  touch Jenkins/dag_fastqs/S8_1.fastq.gz

  printf "\nmultiple peak callers, assemblies and replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv

  printf "\nmultiple peak callers, assemblies and replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv create_trackhub=True

  printf "\nmultiple peak callers, assemblies and replicates - qc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/genrich_macs2.yaml --config samples=../../Jenkins/atac_seq/complex_samples.tsv create_qc_report=True

fi

#if [ $1 = "scatac-seq" ]; then
#
## ATAC-seq workflow (also covers ChIP-seq workflow)
#WF=scATAC_seq
#
#mkdir -p Jenkins/dag_fastqs
#touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
#touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz
#
#printf "\nscatac-seq default\n"
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --dag | dot -Tsvg > graph.svg
#
#
#fi

if [ $1 = "rna-seq" ]; then

# RNA-seq workflow
WF=rna_seq

mkdir -p Jenkins/dag_fastqs
touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz

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
