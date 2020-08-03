# to run these tests locally:
#   python setup.py develop
#   bash ./tests/dag_tests.sh TEST

# check if an argument was passed
if [ -z "$1" ]; then
    echo "No test specified"; exit
fi

# remove the test outputs on exit locally
if [[ $(pwd) != *lib/jenkins* ]]; then
  trap "rm -rf tests/local_test_results" EXIT
fi

CORES=28
function assert_rulecount {
  # check if the DAG (stored with  | tee tests/local_test_results/${1}_dag  ) ran rule $2 exactly $3 times
  val=$(cat tests/local_test_results/${1}_dag | grep -w $2 | cut -f2);

  # check if the rule was found in the DAG at all
  if [ -z "$val" ]; then
    # if specified count is zero, that's OK
    (($3!=0)) && printf "\nrule $2 did not run at all. Exiting.\n\n" && exit 1;
  else
    # else check if the count equals the specified number
    (($val!=$3)) && printf "\nrule $2 ran $val times instead of the expected $3 times. Exiting.\n\n" && exit 1;
  fi
  :
}
mkdir -p tests/local_test_results/fastq
touch tests/local_test_results/fastq/S1_1_R1.fastq.gz
touch tests/local_test_results/fastq/S1_1_R2.fastq.gz
touch tests/local_test_results/fastq/S1_2_R1.fastq.gz
touch tests/local_test_results/fastq/S1_2_R2.fastq.gz
touch tests/local_test_results/fastq/S2_1.fastq.gz
touch tests/local_test_results/fastq/S2_2.fastq.gz
touch tests/local_test_results/fastq/S3_1.fastq.gz
touch tests/local_test_results/fastq/S4_1.fastq.gz
touch tests/local_test_results/fastq/S5_1.fastq.gz
touch tests/local_test_results/fastq/S6_1.fastq.gz
touch tests/local_test_results/fastq/S7_1.fastq.gz
touch tests/local_test_results/fastq/S8_1.fastq.gz

if [ $1 = "alignment" ]; then

  # download workflow
  WF=download_fastq

  printf "\ndownload default\n"
  seq2science run download-fastq -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sra2fastq_PE 1
  assert_rulecount $1 sra2fastq_SE 1

  # alignment workflow
  WF=alignment

  printf "\nalignment default\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_index 1
  assert_rulecount $1 mark_duplicates 1

  printf "\naligners\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:bowtie2,samples:tests/alignment/stranded_sample.tsv,fastq_dir:../tinyfastq} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bowtie2_index 1
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:bwa-mem,samples:tests/alignment/stranded_sample.tsv,fastq_dir:../tinyfastq} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_index 1
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:bwa-mem2,samples:tests/alignment/stranded_sample.tsv,fastq_dir:../tinyfastq} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2_index 1
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:hisat2,samples:tests/alignment/stranded_sample.tsv,fastq_dir:../tinyfastq} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 hisat2_index 1
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:star,samples:tests/alignment/stranded_sample.tsv,fastq_dir:../tinyfastq} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 star_index 1

  printf "\nalignmentsieve\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/alignmentsieve.yaml --snakemakeOptions dryrun=True quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sieve_bam 1
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/alignmentsieve_off.yaml --snakemakeOptions dryrun=True quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sieve_bam 0

  printf "\nsorting\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/samtools_coordinate.yaml --snakemakeOptions dryrun=True quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 samtools_presort 1
  assert_rulecount $1 sambamba_sort 0
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/samtools_queryname.yaml --snakemakeOptions dryrun=True quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 samtools_presort 1
  assert_rulecount $1 sambamba_sort 0
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/sambamba_coordinate.yaml --snakemakeOptions dryrun=True quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sambamba_sort 1
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/sambamba_queryname.yaml --snakemakeOptions dryrun=True quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sambamba_sort 1

  printf "\ncram support\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={cram_no_bam:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam2cram 1

  printf "\ntrackhub\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 1
  assert_rulecount $1 cytoband 1
  assert_rulecount $1 gcPercent 1
  assert_rulecount $1 softmask_track_2 1
  assert_rulecount $1 twobit 1

  printf "\ntrackhub - stranded bams\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={create_trackhub:True,samples:tests/alignment/stranded_sample.tsv,fastq_dir:../tinyfastq} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 2

  printf "\nmultiqc report\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 4

  printf "\nmultiple assemblies\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/assemblies.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_index 2

  printf "\nmultiple assemblies - trackhub\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/assemblies.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 twobit 2

  printf "\nmultiple assemblies - multiqc report\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/assemblies.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 samtools_stats 2

  printf "\nmultiple replicates\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:keep,samples:tests/alignment/replicates.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 bwa_mem 2
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,samples:tests/alignment/replicates.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem 1

  printf "\nmultiple replicates - trackhubs\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,samples:tests/alignment/replicates.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 1

  printf "\nmultiple replicates - multiqc report\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,samples:tests/alignment/replicates.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 8

  printf "\nmultiple assemblies and replicates\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:keep,samples:tests/alignment/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_index 2
  assert_rulecount $1 bwa_mem 4

  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,samples:tests/alignment/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_index 2
  assert_rulecount $1 bwa_mem 2

  printf "\nmultiple assemblies and replicates - trackhub\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,create_trackhub:True,samples:tests/alignment/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 2

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  seq2science run alignment -n --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,create_qc_report:True,samples:tests/alignment/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 samtools_stats 2

  test_ran=1
fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq

  printf "\natac-seq default\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 macs2_callpeak 1

  printf "\npeak callers\n"
  # seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/macs2.yaml --snakemakeOptions dryrun=True quiet=True | tee tests/local_test_results/${1}_dag
  # assert_rulecount $1 macs2_callpeak 1
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich.yaml --snakemakeOptions dryrun=True quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 call_peak_genrich 1

  printf "\ntrackhub\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 1

  printf "\nmultiqc report\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 4

  printf "\nmultiple peak callers\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 macs2_callpeak 1
  assert_rulecount $1 call_peak_genrich 1

  printf "\nmultiple peak callers - trackhub\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 2
  assert_rulecount $1 bedgraphish_to_bedgraph 1

  printf "\nmultiple peak callers - multiqc report\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featureCounts 2
  assert_rulecount $1 coverage_table 2
  assert_rulecount $1 quantile_normalization 2
  assert_rulecount $1 edgeR_normalization 6
  assert_rulecount $1 mean_center 8

  printf "\nmultiple peak callers & multiple assemblies\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/assemblies.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 coverage_table 4

  printf "\nmultiple peak callers & multiple assemblies - trackhub\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/assemblies.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 4

  printf "\nmultiple peak callers & multiple assemblies - multiqc report\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/assemblies.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featureCounts 4

  printf "\nmultiple peak callers & multiple replicates\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/replicates.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem 1

  printf "\nmultiple peak callers & multiple replicates - trackhub\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/replicates.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 2

  printf "\nmultiple peak callers & multiple replicates - multiqc report\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/replicates.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featureCounts 2

  printf "\nmultiple peak callers, assemblies and replicates\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/atac_seq/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem 8
  assert_rulecount $1 coverage_table 4

  printf "\nmultiple peak callers, assemblies and replicates - trackhub\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/atac_seq/complex_samples.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 16

  printf "\nmultiple peak callers, assemblies and replicates - multiqc report\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/atac_seq/complex_samples.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featureCounts 16

  printf "\control and merging of tecnical replicates\n"
  seq2science run atac-seq -n --cores $CORES --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/atac_seq/control.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem 7

  test_ran=1
fi

if [ $1 = "scatac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=scatac_seq

  printf "\nscatac-seq default\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 create_SNAP_object 1

  printf "\ntrackhub\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/samples.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 trackhub 1

  printf "\nqc multiqc report\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/samples.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 4  # twice for sample and twice for trep

  printf "\nmultiple assemblies\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/assemblies.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_index 2
  assert_rulecount $1 create_SNAP_object 2

  printf "\nmultiple assemblies - trackhubs\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/assemblies.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 twobit 2

  printf "\nmultiple assemblies - multiqc\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/assemblies.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 8  # twice for sample and twice for trep

  printf "\nmultiple replicates\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/dag_sample.tsv,technical_replicates:merge} | tee tests/local_test_results/${1}_dag  # nothing to merge
  assert_rulecount $1 merge_replicates 0
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/replicates.tsv,technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 bwa_mem 2
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/replicates.tsv,technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 2
  assert_rulecount $1 bwa_mem 1

  printf "\nmultiple replicates - trackhub\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/replicates.tsv,technical_replicates:merge,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 twobit 1

  printf "\nmultiple replicates - multiqc report\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/replicates.tsv,technical_replicates:merge,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 8  # twice for sample and twice for trep

  printf "\nmultiple assemblies and replicates\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/dag_sample.tsv,technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 bwa_mem 1
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/complex_samples.tsv,technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 bwa_mem 4
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/complex_samples.tsv,technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 3
  assert_rulecount $1 bwa_mem 2

  printf "\nmultiple assemblies and replicates - trackhub\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/complex_samples.tsv,technical_replicates:merge,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 twobit 2

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  seq2science run scatac-seq -n --cores $CORES --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/scatac_seq/complex_samples.tsv,technical_replicates:merge,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 12  # twice for sample and twice for trep

  test_ran=1
fi

if [ $1 = "rna-seq" ]; then

  # RNA-seq workflow
  WF=rna_seq

  printf "\nrna-seq default\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:star} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 htseq_count 1

  printf "\naligners\n"
#  seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:star} | tee tests/local_test_results/${1}_dag  # default
#  assert_rulecount $1 star_align 1
  seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:hisat2} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 hisat2_align 1

  printf "\nquantifiers\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:star,quantifier:salmon} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 1
  # seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:star,quantifier:htseq} | tee tests/local_test_results/${1}_dag  # default
  # assert_rulecount $1 htseq_count 1
  seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=False config={aligner:star,quantifier:featurecounts} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featurecounts 1

  printf "\ndecoy aware salmon index\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/salmon_config.yaml --snakemakeOptions dryrun=True quiet=True config={samples:tests/alignment/dag_sample.tsv,fastq_dir:tests/local_test_results/fastq} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 decoy_transcripts 1

  printf "\ntrackhub\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:star,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 0
  assert_rulecount $1 bam_bigwig 1
  seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:star,quantifier:salmon,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 1
  assert_rulecount $1 bam_bigwig 1

  printf "\nmultiqc report\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions dryrun=True quiet=True config={aligner:star,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 4

  printf "\ndifferential expression analysis\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 htseq_count 10
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={quantifier:salmon,technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 10

  printf "\nmultiple assemblies with DEA\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:keep,samples:tests/rna_seq/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 star_index 2
  assert_rulecount $1 star_align 10

  printf "\nmultiple assemblies with DEA - trackhubs\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:keep,samples:tests/rna_seq/complex_samples.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 10

  printf "\nmultiple assemblies with DEA - multiqc\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:keep,samples:tests/rna_seq/complex_samples.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 24
  assert_rulecount $1 multiqc 2

  printf "\nmultiple replicates with DEA \n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 htseq_count 10
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 htseq_count 8

  printf "\nmultiple replicates with DEA - trackhubs\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 8

  printf "\nmultiple replicates with DEA - multiqc\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 24
  assert_rulecount $1 htseq_count 8

  printf "\nmultiple assemblies and replicates with DEA \n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:keep,samples:tests/rna_seq/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 htseq_count 10
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,samples:tests/rna_seq/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 htseq_count 8
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={quantifier:salmon,technical_replicates:merge,samples:tests/rna_seq/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 8

  printf "\nmultiple assemblies and replicates with DEA - trackhub\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={quantifier:salmon,technical_replicates:merge,samples:tests/rna_seq/complex_samples.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 8
  assert_rulecount $1 star_align 8
  assert_rulecount $1 bam_bigwig 8

  printf "\nmultiple assemblies and replicates with DEA - multiqc report\n"
  seq2science run rna-seq -n --cores $CORES --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions dryrun=True quiet=True config={technical_replicates:merge,samples:tests/rna_seq/complex_samples.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc  24

  test_ran=1
fi

# check if any test has run
if [ -z "$test_ran" ]; then
  printf "\nunrecognized input: ${1}\n"; exit
fi
