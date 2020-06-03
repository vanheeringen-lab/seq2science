#!/usr/bin/env bash

# to run these tests locally:
# 1)   chmod u+x ./tests/dag_tests.sh
# 2)   ./tests/dag_tests.sh TEST

if [ -z "$1" ]
  then
    echo "No test specified"
    exit
fi

CORES=48
trap "rm -rf Jenkins_results" EXIT  # remove the test outputs on exit
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
mkdir -p Jenkins_results/fastq
touch Jenkins_results/fastq/S1_1_R1.fastq.gz
touch Jenkins_results/fastq/S1_1_R2.fastq.gz
touch Jenkins_results/fastq/S1_2_R1.fastq.gz
touch Jenkins_results/fastq/S1_2_R2.fastq.gz
touch Jenkins_results/fastq/S2_1.fastq.gz
touch Jenkins_results/fastq/S2_2.fastq.gz
touch Jenkins_results/fastq/S3_1.fastq.gz
touch Jenkins_results/fastq/S4_1.fastq.gz
touch Jenkins_results/fastq/S5_1.fastq.gz
touch Jenkins_results/fastq/S6_1.fastq.gz
touch Jenkins_results/fastq/S7_1.fastq.gz
touch Jenkins_results/fastq/S8_1.fastq.gz

if [ $1 = "alignment" ]; then

  # download workflow
  WF=download_fastq

  printf "\ndownload default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml | tee Jenkins_results/val
  assert_rulecount sra2fastq_PE 1
  assert_rulecount sra2fastq_SE 1

  # alignment workflow
  WF=alignment

  printf "\nalignment default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml | tee Jenkins_results/val
  assert_rulecount bwa_index 1
  assert_rulecount mark_duplicates 1

  printf "\ntrackhub - stranded bams\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/stranded_sample.tsv create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bam_stranded_bigwig 1

  printf "\naligners\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml  --config aligner=bowtie2 | tee Jenkins_results/val
  assert_rulecount bowtie2_index 1
#  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml  --config aligner=bwa | tee Jenkins_results/val  # default
#  assert_rulecount bwa_index 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml  --config aligner=hisat2 | tee Jenkins_results/val
  assert_rulecount hisat2_index 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml  --config aligner=star | tee Jenkins_results/val
  assert_rulecount star_index 1

  printf "\nalignmentsieve\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/alignmentsieve.yaml | tee Jenkins_results/val
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/alignmentsieve_off.yaml | tee Jenkins_results/val

  printf "\nsorting\n"
  # snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/samtools_coordinate.yaml  # default
  # assert_rulecount samtools_presort 1
  # assert_rulecount sambamba_sort 0
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/samtools_queryname.yaml | tee Jenkins_results/val
  assert_rulecount samtools_sort 1
  assert_rulecount sambamba_sort 0
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/sambamba_coordinate.yaml | tee Jenkins_results/val
  assert_rulecount sambamba_sort 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/sambamba_queryname.yaml | tee Jenkins_results/val
  assert_rulecount sambamba_sort 1

  printf "\ncram support\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config cram_no_bam=True | tee Jenkins_results/val
  assert_rulecount bam2cram 1

  printf "\ntrackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bam_bigwig 1
  assert_rulecount cytoband 1
  assert_rulecount gcPercent 1
  assert_rulecount softmask_track_2 1
  assert_rulecount twobit 1

  printf "\nmultiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 4

  printf "\nmultiple assemblies\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/assemblies.tsv | tee Jenkins_results/val
  assert_rulecount bwa_index 2

  printf "\nmultiple assemblies - trackhubs\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/assemblies.tsv create_trackhub=True | tee Jenkins_results/val
  assert_rulecount twobit 2

  printf "\nmultiple assemblies - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/assemblies.tsv create_qc_report=True | tee Jenkins_results/val
  assert_rulecount samtools_stats 2

  printf "\nmultiple replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config technical_replicates=merge | tee Jenkins_results/val  # nothing to merge
  assert_rulecount merge_replicates 0
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/replicates.tsv technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount merge_replicates 0
  assert_rulecount bwa_mem 2
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/replicates.tsv technical_replicates=merge | tee Jenkins_results/val
  assert_rulecount bwa_mem 1

  printf "\nmultiple replicates - trackhubs\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/replicates.tsv technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bam_bigwig 1

  printf "\nmultiple replicates - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/replicates.tsv technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 8

  printf "\nmultiple assemblies and replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/complex_samples.tsv technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount bwa_index 2
  assert_rulecount bwa_mem 4

  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/complex_samples.tsv technical_replicates=merge | tee Jenkins_results/val
  assert_rulecount bwa_index 2
  assert_rulecount bwa_mem 2

  printf "\nmultiple assemblies and replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/complex_samples.tsv technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bam_bigwig 2

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/default_config.yaml --config samples=../../tests/alignment/complex_samples.tsv technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
  assert_rulecount samtools_stats 2

fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq

  printf "\natac-seq default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml | tee Jenkins_results/val
  assert_rulecount macs2_callpeak 1

  printf "\npeak callers\n"
  # snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/macs2.yaml  | tee Jenkins_results/val  # default
  #  assert_rulecount macs2_callpeak 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich.yaml | tee Jenkins_results/val
  assert_rulecount call_peak_genrich 1

  printf "\ntrackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bedgraph_bigwig 1

  printf "\nmultiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 4

  printf "\nmultiple peak callers\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml | tee Jenkins_results/val
  assert_rulecount macs2_callpeak 1
  assert_rulecount call_peak_genrich 1

  printf "\nmultiple peak callers - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bedgraph_bigwig 2
  assert_rulecount bedgraphish_to_bedgraph 1

  printf "\nmultiple peak callers - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config create_qc_report=True | tee Jenkins_results/val
  assert_rulecount featureCounts 2

  printf "\nmultiple peak callers & multiple assemblies\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/alignment/assemblies.tsv | tee Jenkins_results/val
  assert_rulecount coverage_table 4

  printf "\nmultiple peak callers & multiple assemblies - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/alignment/assemblies.tsv create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bedgraph_bigwig 4

  printf "\nmultiple peak callers & multiple assemblies - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/alignment/assemblies.tsv create_qc_report=True | tee Jenkins_results/val
  assert_rulecount featureCounts 4

  printf "\nmultiple peak callers & multiple replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/alignment/replicates.tsv | tee Jenkins_results/val
  assert_rulecount bwa_mem 1

  printf "\nmultiple peak callers & multiple replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/alignment/replicates.tsv create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bedgraph_bigwig 2

  printf "\nmultiple peak callers & multiple replicates - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/alignment/replicates.tsv create_qc_report=True | tee Jenkins_results/val
  assert_rulecount featureCounts 2

  printf "\nmultiple peak callers, assemblies and replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/atac_seq/complex_samples.tsv | tee Jenkins_results/val
  assert_rulecount bwa_mem 8
  assert_rulecount coverage_table 4

  printf "\nmultiple peak callers, assemblies and replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/atac_seq/complex_samples.tsv create_trackhub=True | tee Jenkins_results/val
  assert_rulecount bedgraph_bigwig 16

  printf "\nmultiple peak callers, assemblies and replicates - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/genrich_macs2.yaml --config samples=../../tests/atac_seq/complex_samples.tsv create_qc_report=True | tee Jenkins_results/val
  assert_rulecount featureCounts 16

fi

if [ $1 = "scatac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=scATAC_seq

  printf "\nscatac-seq default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml | tee Jenkins_results/val
  assert_rulecount create_SNAP_object 1

  printf "\ntrackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config create_trackhub=True | tee Jenkins_results/val
  # TODO: scATAC-seq does not create a trackhub
  #assert_rulecount trackhub 1

  printf "\nqc multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 2

  printf "\nmultiple assemblies\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/assemblies.tsv | tee Jenkins_results/val
  assert_rulecount bwa_index 2
  assert_rulecount create_SNAP_object 2

  printf "\nmultiple assemblies - trackhubs\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/assemblies.tsv create_trackhub=True | tee Jenkins_results/val
  # TODO: scATAC-seq does not create a trackhub
  #assert_rulecount bam_bigwig 2
  #assert_rulecount twobit 2

  printf "\nmultiple assemblies - multiqc\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/assemblies.tsv create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 4

  printf "\nmultiple replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config technical_replicates=merge | tee Jenkins_results/val  # nothing to merge
  assert_rulecount merge_replicates 0
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/replicates.tsv technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount merge_replicates 0
  assert_rulecount bwa_mem 2
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/replicates.tsv technical_replicates=merge | tee Jenkins_results/val
  assert_rulecount merge_replicates 2
  assert_rulecount bwa_mem 1

  printf "\nmultiple replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/replicates.tsv technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
  # TODO: scATAC-seq does not create a trackhub
  #assert_rulecount bam_bigwig 1

  printf "\nmultiple replicates - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/replicates.tsv technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
  # different number from other workflows
  assert_rulecount fastqc 2

  printf "\nmultiple assemblies and replicates\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/complex_samples.tsv technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount merge_replicates 0
  assert_rulecount bwa_mem 4
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/complex_samples.tsv technical_replicates=merge | tee Jenkins_results/val
  assert_rulecount merge_replicates 3
  assert_rulecount bwa_mem 2

  printf "\nmultiple assemblies and replicates - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/complex_samples.tsv technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
  # TODO: scATAC-seq does not create a trackhub
  #assert_rulecount bam_bigwig 2

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config samples=../../tests/alignment/complex_samples.tsv technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
  # different number from other workflows
  assert_rulecount fastqc 4

fi

if [ $1 = "rna-seq" ]; then

  # RNA-seq workflow
  WF=rna_seq

  printf "\nrna-seq default\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config quantifier=star | tee Jenkins_results/val
  assert_rulecount star_quant 1

  printf "\nquantifiers\n"
  # snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config quantifier=star | tee Jenkins_results/val  # default
  # assert_rulecount star_quant 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config quantifier=salmon | tee Jenkins_results/val
  assert_rulecount salmon_quant 1

  printf "\ndecoy aware salmon index\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config quantifier=salmon decoy_aware_index=True | tee Jenkins_results/val
  # TODO: bug: decoy not used!
  #assert_rulecount decoy_transcripts 1

  printf "\ntrackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config quantifier=star create_trackhub=True | tee Jenkins_results/val
  assert_rulecount salmon_quant 0
  assert_rulecount star_quant 0
  assert_rulecount star_align 1
  assert_rulecount bam_bigwig 1
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config quantifier=salmon create_trackhub=True | tee Jenkins_results/val
  assert_rulecount salmon_quant 1
  assert_rulecount star_quant 0
  assert_rulecount star_align 1
  assert_rulecount bam_bigwig 1

  printf "\nmultiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/alignment/default_config.yaml --config quantifier=star create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 4

  printf "\ndifferential expression analysis\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config quantifier=star technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount star_quant 10
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config quantifier=salmon technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount salmon_quant 10

  printf "\nmultiple assemblies with DEA\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config samples=../../tests/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount star_index 2
  # TODO: bug: quantifier runs 2x too many times (2 assemblies)
  #assert_rulecount star_quant 10

  printf "\nmultiple assemblies with DEA - trackhubs\n"
  # TODO: error!
  #snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config samples=../../tests/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep create_trackhub=True | tee Jenkins_results/val
  #assert_rulecount bam_bigwig 20

  printf "\nmultiple assemblies with DEA - multiqc\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config samples=../../tests/rna_seq/complex_samples.tsv quantifier=star technical_replicates=keep create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 24
  assert_rulecount multiqc 2

  printf "\nmultiple replicates with DEA \n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config quantifier=star technical_replicates=keep | tee Jenkins_results/val
  assert_rulecount merge_replicates 0
  assert_rulecount star_quant 10
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config quantifier=star technical_replicates=merge | tee Jenkins_results/val
  assert_rulecount star_quant 8

  printf "\nmultiple replicates with DEA - trackhubs\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config quantifier=star technical_replicates=merge create_trackhub=True | tee Jenkins_results/val
  # TODO: error!
  #assert_rulecount bam_bigwig 8

  printf "\nmultiple replicates with DEA - multiqc\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config quantifier=star technical_replicates=merge create_qc_report=True | tee Jenkins_results/val
  assert_rulecount fastqc 24
  assert_rulecount star_quant 8

  printf "\nmultiple assemblies and replicates with DEA \n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config samples=../../tests/rna_seq/complex_samples.tsv technical_replicates=keep quantifier=star | tee Jenkins_results/val
  assert_rulecount merge_replicates 0
  # TODO: bug: quantifier runs 2x too many times (2 assemblies)
  #assert_rulecount star_quant 10
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config samples=../../tests/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=star | tee Jenkins_results/val
  # TODO: bug: quantifier runs 2x too many times (2 assemblies)
  #assert_rulecount star_quant 8
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config samples=../../tests/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=salmon | tee Jenkins_results/val
  # TODO: bug: quantifier runs 2x too many times (2 assemblies)
  #assert_rulecount salmon_quant 8

  printf "\nmultiple assemblies and replicates with DEA - trackhub\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config samples=../../tests/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=salmon create_trackhub=True | tee Jenkins_results/val
  # TODO: bug: quantifier runs 16 (2x8) times, aligner runs 8 times.
  #assert_rulecount salmon_quant 8
  assert_rulecount star_align 8
  assert_rulecount bam_bigwig 8

  printf "\nmultiple assemblies and replicates with DEA - multiqc report\n"
  snakemake -n -j $CORES --quiet -s workflows/$WF/Snakefile --directory workflows/$WF --configfile tests/$WF/rna_seq_config.yaml --config samples=../../tests/rna_seq/complex_samples.tsv technical_replicates=merge quantifier=star create_qc_report=True | tee Jenkins_results/val
  # TODO: bug: quantifier runs 16 (2x8) times, aligner runs 8 times.
  assert_rulecount fastqc  24

fi
