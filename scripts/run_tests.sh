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
# TODO: trap "rm -rf Jenkins_results;" EXIT  # remove the test outputs on exit
set -e  # Exit immediately if a command exits with a non-zero status.
#function assert_rulecount {
#  # check if the DAG (stored with  | tee Jenkins_results/val  ) ran rule $1 exactly $2 times
#  val=$(cat Jenkins_results/val | grep -w $1 | cut -f2);
#
#  # check if the rule was found in the DAG at all
#  if [ -z "$val" ]; then
#    # if specified count is zero, that's OK
#    (($2!=0)) && printf "\nrule $1 did not run at all. Exiting.\n\n" && exit 1;
#  else
#    # else check if the count equals the specified number
#    (($val!=$2)) && printf "\nrule $1 ran $val times instead of the expected $2 times. Exiting.\n\n" && exit 1;
#  fi
#  :
#}
#mkdir -p Jenkins_results

if [ $1 = "download" ]; then

  WF=download_fastq

  # test basic downloading 1 PE and 1 SE
  printf "\ndownload SE and PE fastqs\n\n"
  snakemake --use-conda -j $CORES -s workflows/$WF/Snakefile --directory workflows/$WF \
  --configfile Jenkins/$WF/default_config.yaml \
  --config samples=../../Jenkins/download_fastq/remote_samples.tsv

  WF=alignment

  # test genome & annotation downloading
  printf "\ndownload genome & annotation\n\n"
  snakemake --use-conda -j $CORES -s workflows/$WF/Snakefile --directory workflows/$WF \
  --configfile Jenkins/$WF/remote_genome_n_sample.yaml \
  -O id2sra sra2fastq_SE sra2fastq_PE fastqc trim_galore_PE get_annotation complement_blacklist setup_blacklist

fi

if [ $1 = "prep_align" ]; then

  snakemake -s workflows/alignment/Snakefile --directory workflows/alignment \
  --use-conda -j $CORES --configfile Jenkins/alignment/default_config.yaml \
  --config samples=../../Jenkins/alignment/remote_genome_n_sample.tsv aligner=bowtie2 \
  --create-envs-only

  snakemake -s workflows/alignment/Snakefile --directory workflows/alignment \
  --use-conda -j $CORES --configfile Jenkins/alignment/default_config.yaml \
  --config samples=../../Jenkins/alignment/remote_genome_n_sample.tsv aligner=bwa \
  --create-envs-only

  snakemake -s workflows/alignment/Snakefile --directory workflows/alignment \
  --use-conda -j $CORES --configfile Jenkins/alignment/default_config.yaml \
  --config samples=../../Jenkins/alignment/remote_genome_n_sample.tsv aligner=hisat2 \
  --create-envs-only

  snakemake -s workflows/alignment/Snakefile --directory workflows/alignment \
  --use-conda -j $CORES --configfile Jenkins/alignment/default_config.yaml \
  --config samples=../../Jenkins/alignment/remote_genome_n_sample.tsv aligner=star \
  --create-envs-only

fi

if [ $1 = "bowtie2" ]; then

  ALIGNER=bowtie2
  WF=alignment
  let "c = $CORES / 4"
  let "a = $c - 2"
  let "s = 2"

  snakemake -s workflows/$WF/Snakefile --directory workflows/$WF \
  --use-conda --nolock --notemp \
  --configfile \
      Jenkins/$WF/default_config.yaml \
  --config \
      aligner=$ALIGNER \
      samples=../../Jenkins/alignment/local_sample.tsv \
      fastq_dir=../../Jenkins/tinyfastq \
      genome_dir=../../Jenkins \
      result_dir=../../Jenkins_results/$ALIGNER \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

fi

if [ $1 = "bwa" ]; then

  ALIGNER=bwa
  WF=alignment
  let "c = $CORES / 4"
  let "a = $c - 2"
  let "s = 2"

  snakemake -s workflows/$WF/Snakefile --directory workflows/$WF \
  --use-conda --nolock --notemp \
  --configfile \
      Jenkins/$WF/default_config.yaml \
  --config \
      aligner=$ALIGNER \
      samples=../../Jenkins/alignment/local_sample.tsv \
      fastq_dir=../../Jenkins/tinyfastq \
      genome_dir=../../Jenkins \
      result_dir=../../Jenkins_results/$ALIGNER \
  -j $c --set-threads bwa_mem=$a samtools_presort=$s

fi

if [ $1 = "hisat2" ]; then

  ALIGNER=hisat2
  WF=alignment
  let "c = $CORES / 4"
  let "a = $c - 2"
  let "s = 2"

  snakemake -s workflows/$WF/Snakefile --directory workflows/$WF \
  --use-conda --nolock --notemp \
  --configfile \
      Jenkins/$WF/default_config.yaml \
  --config \
      aligner=$ALIGNER \
      samples=../../Jenkins/alignment/local_sample.tsv \
      fastq_dir=../../Jenkins/tinyfastq \
      genome_dir=../../Jenkins \
      result_dir=../../Jenkins_results/$ALIGNER \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

fi

if [ $1 = "star" ]; then

  ALIGNER=star
  WF=alignment
  let "c = $CORES / 4"
  let "a = $c - 2"
  let "s = 2"

  snakemake -s workflows/$WF/Snakefile --directory workflows/$WF \
  --use-conda --nolock --notemp \
  --configfile \
      Jenkins/$WF/default_config.yaml \
  --config \
      aligner=$ALIGNER \
      samples=../../Jenkins/alignment/local_sample.tsv \
      fastq_dir=../../Jenkins/tinyfastq \
      genome_dir=../../Jenkins \
      result_dir=../../Jenkins_results/$ALIGNER \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq
  ALIGNER=bowtie2

#  printf "\natac-seq defaul\n"
#  snakemake -s workflows/$WF/Snakefile --directory workflows/$WF \
#  --use-conda -j $CORES \
#  --configfile \
#      Jenkins/alignment/remote_genome_n_sample.yaml \
#  --config \
#      aligner=bowtie2

  printf "\natac-seq with multiqc report\n"
  # creating both at once because the whole workflow needs to run again, and the remote sample is slow
  snakemake -s workflows/$WF/Snakefile --directory workflows/$WF \
  --use-conda -j $CORES \
  --configfile \
      Jenkins/alignment/remote_genome_n_sample.yaml \
  --config \
      aligner=bowtie2 \
      create_qc_report=True

  printf "\ntrackhub\n"
  snakemake -s workflows/$WF/Snakefile --directory workflows/$WF \
  --use-conda -j $CORES \
  --configfile \
      Jenkins/alignment/remote_genome_n_sample.yaml \
  --config \
      aligner=bowtie2 \
      create_trackhub=True

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
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config create_qc_report=True | tee Jenkins_results/val
assert_rulecount fastqc 2

printf "\nmultiple assemblies\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv | tee Jenkins_results/val
assert_rulecount bwa_index 2
assert_rulecount create_SNAP_object 2

printf "\nmultiple assemblies - trackhubs\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_trackhub=True | tee Jenkins_results/val
# TODO: scATAC-seq does not create a trackhub
#assert_rulecount bam_bigwig 2
#assert_rulecount twobit 2

printf "\nmultiple assemblies - multiqc\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/assemblies.tsv create_qc_report=True | tee Jenkins_results/val
assert_rulecount fastqc 4

printf "\nmultiple replicates\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config technical_replicates=merge | tee Jenkins_results/val  # nothing to merge
assert_rulecount merge_replicates 0
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=keep | tee Jenkins_results/val
assert_rulecount merge_replicates 0
assert_rulecount bwa_mem 2
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=merge | tee Jenkins_results/val
assert_rulecount merge_replicates 2
assert_rulecount bwa_mem 1

printf "\nmultiple replicates - trackhub\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
# TODO: scATAC-seq does not create a trackhub
#assert_rulecount bam_bigwig 1

printf "\nmultiple replicates - multiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/replicates.tsv technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
# different number from other workflows
assert_rulecount fastqc 2

printf "\nmultiple assemblies and replicates\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=keep | tee Jenkins_results/val
assert_rulecount merge_replicates 0
assert_rulecount bwa_mem 4
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge | tee Jenkins_results/val
assert_rulecount merge_replicates 3
assert_rulecount bwa_mem 2

printf "\nmultiple assemblies and replicates - trackhub\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
# TODO: scATAC-seq does not create a trackhub
#assert_rulecount bam_bigwig 2

printf "\nmultiple assemblies and replicates - multiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config samples=../../Jenkins/alignment/complex_samples.tsv technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
# different number from other workflows
assert_rulecount fastqc 4

fi

if [ $1 = "rna-seq" ]; then

# RNA-seq workflow
WF=rna_seq

#mkdir -p Jenkins/dag_fastqs
#touch Jenkins/dag_fastqs/S1_1_R1.fastq.gz
#touch Jenkins/dag_fastqs/S1_1_R2.fastq.gz

printf "\nrna-seq default\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=star | tee Jenkins_results/val
assert_rulecount star_quant 1

printf "\nquantifiers\n"
# snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=star | tee Jenkins_results/val  # default
# assert_rulecount star_quant 1
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=salmon | tee Jenkins_results/val
assert_rulecount salmon_quant 1

printf "\ndecoy aware salmon index\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=salmon decoy_aware_index=True | tee Jenkins_results/val
# TODO: bug: decoy not used!
#assert_rulecount decoy_transcripts 1

printf "\ntrackhub\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=star create_trackhub=True | tee Jenkins_results/val
assert_rulecount salmon_quant 0
assert_rulecount star_quant 0
assert_rulecount star_align 1
assert_rulecount bam_bigwig 1
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=salmon create_trackhub=True | tee Jenkins_results/val
assert_rulecount salmon_quant 1
assert_rulecount star_quant 0
assert_rulecount star_align 1
assert_rulecount bam_bigwig 1

printf "\nmultiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/alignment/default_config.yaml --config quantifier=star create_qc_report=True | tee Jenkins_results/val
assert_rulecount fastqc 4

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
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config quantifier=star technical_replicates=keep | tee Jenkins_results/val
assert_rulecount star_quant 10
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config quantifier=salmon technical_replicates=keep | tee Jenkins_results/val
assert_rulecount salmon_quant 10

printf "\nmultiple assemblies with DEA\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep | tee Jenkins_results/val
assert_rulecount star_index 2
# TODO: bug: quantifier runs 2x too many times (2 assemblies)
#assert_rulecount star_quant 10

printf "\nmultiple assemblies with DEA - trackhubs\n"
# TODO: error!
#snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep create_trackhub=True | tee Jenkins_results/val
#assert_rulecount bam_bigwig 20

printf "\nmultiple assemblies with DEA - multiqc\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep create_qc_report=True | tee Jenkins_results/val
assert_rulecount fastqc 24
assert_rulecount multiqc 2

printf "\nmultiple replicates with DEA \n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config quantifier=star technical_replicates=keep | tee Jenkins_results/val
assert_rulecount merge_replicates 0
assert_rulecount star_quant 10
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config quantifier=star technical_replicates=merge | tee Jenkins_results/val
assert_rulecount star_quant 8

printf "\nmultiple replicates with DEA - trackhubs\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config quantifier=star technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
# TODO: error!
#assert_rulecount bam_bigwig 8

printf "\nmultiple replicates with DEA - multiqc\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config quantifier=star technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
assert_rulecount fastqc 24
assert_rulecount star_quant 8

printf "\nmultiple assemblies and replicates with DEA \n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=keep quantifier=star | tee Jenkins_results/val
assert_rulecount merge_replicates 0
# TODO: bug: quantifier runs 2x too many times (2 assemblies)
#assert_rulecount star_quant 10
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=star | tee Jenkins_results/val
# TODO: bug: quantifier runs 2x too many times (2 assemblies)
#assert_rulecount star_quant 8
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=salmon | tee Jenkins_results/val
# TODO: bug: quantifier runs 2x too many times (2 assemblies)
#assert_rulecount salmon_quant 8

printf "\nmultiple assemblies and replicates with DEA - trackhub\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=salmon create_trackhub=True | tee Jenkins_results/val
# TODO: bug: quantifier runs 16 (2x8) times, aligner runs 8 times.
#assert_rulecount salmon_quant 8
assert_rulecount star_align 8
assert_rulecount bam_bigwig 8

printf "\nmultiple assemblies and replicates with DEA - multiqc report\n"
snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile Jenkins/$WF/deseq2.yaml --config samples=../../Jenkins/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=star create_qc_report=True | tee Jenkins_results/val
# TODO: bug: quantifier runs 16 (2x8) times, aligner runs 8 times.
assert_rulecount fastqc  24

fi
