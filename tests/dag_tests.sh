# to run these tests locally:
#   mamba install toml --yes
#   python setup.py develop
#   bash ./tests/dag_tests.sh TEST

# for lazy asses:
# bash ./tests/dag_tests.sh alignment && bash ./tests/dag_tests.sh atac-seq && bash ./tests/dag_tests.sh scatac-seq && bash ./tests/dag_tests.sh rna-seq && bash ./tests/dag_tests.sh scrna-seq

# check if an argument was passed
if [ -z "$1" ]; then
    echo "No test specified"; exit
fi

# remove the test outputs on exit locally
if [[ $(pwd) != *lib/jenkins* ]]; then
  trap "rm -rf tests/local_test_results" EXIT
fi

function assert_rulecount {
  # check if the DAG (stored with  | tee tests/local_test_results/${1}_dag  ) ran rule $2 exactly $3 times
  val=$(cat tests/local_test_results/${1}_dag | grep -wE $2$ | cut -f2);
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
touch tests/local_test_results/fastq/S1_2_R1_100.fastq.gz  # Illumina '_100' suffix -> must be renamed
touch tests/local_test_results/fastq/S1_2_R2_100.fastq.gz  # Illumina '_100' suffix -> must be renamed
touch tests/local_test_results/fastq/S2_1.fastq.gz
touch tests/local_test_results/fastq/S2_2.fastq.gz
touch tests/local_test_results/fastq/S3_1_100.fastq.gz  # Illumina '_100' suffix -> must be renamed
touch tests/local_test_results/fastq/S4_1.fastq.gz
touch tests/local_test_results/fastq/S5_1.fastq.gz
touch tests/local_test_results/fastq/S6_1.fastq.gz
touch tests/local_test_results/fastq/S7_1.fastq.gz
touch tests/local_test_results/fastq/S8_1.fastq.gz
touch tests/local_test_results/ERCC92.fa
touch tests/local_test_results/ERCC92.gtf
touch tests/local_test_results/barcodes.txt

for assembly in assembly1 assembly2; do
  mkdir -p tests/local_test_results/${assembly}
  touch tests/local_test_results/${assembly}/${assembly}.fa
  touch tests/local_test_results/${assembly}/${assembly}.annotation.gtf
  touch tests/local_test_results/${assembly}/${assembly}.annotation.bed
done


if [ $1 = "alignment" ]; then

  # download workflow
  WF=download_fastq

  printf "\ndownload default\n"
  seq2science run download-fastq -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 runs2sample 3
  assert_rulecount $1 'ena2fastq_SE|sra2fastq_SE' 5
  assert_rulecount $1 'ena2fastq_PE|sra2fastq_PE' 1

  # alignment workflow
  WF=alignment

  printf "\nalignment default\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2_index 1
  assert_rulecount $1 mark_duplicates 1

  printf "\ntrimmers\n"
  printf "\tfastp\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={trimmer:fastp,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastp_PE 1
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/replicates.tsv,trimmer:fastp,create_qc_report:True,technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 rename_sample 2
  assert_rulecount $1 fastp_PE 2
  assert_rulecount $1 fastp_qc_PE 1

  printf "\n\ttrimgalore\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={trimmer:trimgalore,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 trimgalore_PE 1
  assert_rulecount $1 fastqc 4
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/replicates.tsv,trimmer:trimgalore,create_qc_report:True,technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 trimgalore_PE 2
  assert_rulecount $1 fastqc 8

  printf "\naligners\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={aligner:bowtie2} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bowtie2_index 1
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={aligner:bwa-mem} | tee tests/local_test_results/${1}_dag  # default
  assert_rulecount $1 bwa_index 1
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={aligner:bwa-mem2} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2_index 1
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={aligner:hisat2} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 hisat2_index 1
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={aligner:minimap2} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 minimap2_index 1
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={aligner:star} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 star_index 1

  printf "\nalignmentsieve\n"
  seq2science run alignment -n --configfile tests/$WF/alignmentsieve.yaml --snakemakeOptions quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sieve_bam 1
  seq2science run alignment -n --configfile tests/$WF/alignmentsieve_off.yaml --snakemakeOptions quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sieve_bam 0

  printf "\nsorting\n"
  seq2science run alignment -n --configfile tests/$WF/samtools_coordinate.yaml --snakemakeOptions quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 samtools_presort 1
  assert_rulecount $1 sambamba_sort 0
  seq2science run alignment -n --configfile tests/$WF/samtools_queryname.yaml --snakemakeOptions quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 samtools_presort 1
  assert_rulecount $1 sambamba_sort 0
  seq2science run alignment -n --configfile tests/$WF/sambamba_coordinate.yaml --snakemakeOptions quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sambamba_sort 1
  seq2science run alignment -n --configfile tests/$WF/sambamba_queryname.yaml --snakemakeOptions quiet=True| tee tests/local_test_results/${1}_dag
  assert_rulecount $1 sambamba_sort 1

  printf "\ncram support\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={cram_no_bam:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam2cram 1

  printf "\ntrackhub\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 1
  assert_rulecount $1 cytoband 1
  assert_rulecount $1 gcPercent 1
  assert_rulecount $1 softmask_track_2 1
  assert_rulecount $1 twobit 1

  printf "\nmultiqc report\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 4

  printf "\nmultiple assemblies\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/assemblies.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2_index 2

  printf "\nmultiple assemblies - trackhub\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/assemblies.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 twobit 2

  printf "\nmultiple assemblies - multiqc report\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/assemblies.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 samtools_stats 2

  printf "\nmultiple replicates\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:keep,samples:tests/alignment/replicates.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 bwa_mem2 2
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,samples:tests/alignment/replicates.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2 1
  assert_rulecount $1 fastp_PE 2

  printf "\nmultiple replicates - trackhubs\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,samples:tests/alignment/replicates.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 1

  printf "\nmultiple replicates - multiqc report\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,samples:tests/alignment/replicates.tsv,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 8

  printf "\nmultiple assemblies and replicates\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:keep,samples:tests/alignment/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2_index 2
  assert_rulecount $1 bwa_mem2 4

  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,samples:tests/alignment/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2_index 2
  assert_rulecount $1 bwa_mem2 2

  printf "\nmultiple assemblies and replicates - trackhub\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,create_trackhub:True,samples:tests/alignment/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bam_bigwig 2

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  seq2science run alignment -n --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,create_qc_report:True,samples:tests/alignment/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 samtools_stats 2

  test_ran=1
fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq

  printf "\natac-seq default\n"
  seq2science run atac-seq -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 macs2_callpeak 1

  printf "\npeak callers\n"
  # seq2science run atac-seq -n --configfile tests/$WF/macs2.yaml --snakemakeOptions quiet=True | tee tests/local_test_results/${1}_dag  # default
  # assert_rulecount $1 macs2_callpeak 1
  seq2science run atac-seq -n --configfile tests/$WF/genrich.yaml --snakemakeOptions quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 call_peak_genrich 1

  printf "\ndifferential peak analysis\n"
  seq2science run atac-seq --skip-rerun -n --configfile tests/rna_seq/deseq2_config.yaml --snakemakeOptions quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 coverage_table 1
  assert_rulecount $1 deseq2 1

  printf "\ntrackhub\n"
  seq2science run atac-seq -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 1

  printf "\nmultiqc report\n"
  seq2science run atac-seq -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 4

  printf "\natac-seq template length filter\n"
  seq2science run atac-seq -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 mark_duplicates 1

  seq2science run atac-seq -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={create_qc_report:True,min_template_length:150} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 mark_duplicates 2

  seq2science run atac-seq -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={create_qc_report:True,max_template_length:150} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 mark_duplicates 2

  printf "\ncustom assembly\n"
  seq2science run atac-seq -n --configfile tests/$WF/macs2.yaml --snakemakeOptions quiet=True config={custom_genome_extension:tests/local_test_results/ERCC92.fa} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 extend_genome 1
  assert_rulecount $1 bwa_mem2_index 1
  assert_rulecount $1 bwa_mem2 1
  seq2science run atac-seq -n --configfile tests/$WF/genrich.yaml --snakemakeOptions quiet=True config={custom_genome_extension:tests/local_test_results/ERCC92.fa} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 extend_genome 1
  assert_rulecount $1 bwa_mem2_index 1
  assert_rulecount $1 bwa_mem2 1

  printf "\nmultiple peak callers\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 macs2_callpeak 1
  assert_rulecount $1 call_peak_genrich 1

  printf "\nmultiple peak callers - trackhub\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 2
  assert_rulecount $1 bedgraphish_to_bedgraph 1

  printf "\nmultiple peak callers - multiqc report\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featureCounts 2
  assert_rulecount $1 coverage_table 2
  assert_rulecount $1 quantile_normalization 2
  assert_rulecount $1 edgeR_normalization 6
  assert_rulecount $1 mean_center 8

  printf "\nmultiple peak callers & multiple assemblies\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/assemblies.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 coverage_table 4

  printf "\nmultiple peak callers & multiple assemblies - trackhub\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/assemblies.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 4

  printf "\nmultiple peak callers & multiple assemblies - multiqc report\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/assemblies.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featureCounts 4

  printf "\nmultiple peak callers & multiple replicates\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/replicates.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2 1

  printf "\nmultiple peak callers & multiple replicates - trackhub\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/replicates.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 2

  printf "\nmultiple peak callers & multiple replicates - multiqc report\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/replicates.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featureCounts 2

  printf "\nmultiple peak callers, assemblies and replicates\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/atac_seq/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 rename_sample 3
  assert_rulecount $1 bwa_mem2 8
  assert_rulecount $1 coverage_table 4

  printf "\nmultiple peak callers, assemblies and replicates - trackhub\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/atac_seq/complex_samples.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bedgraph_bigwig 16

  printf "\nmultiple peak callers, assemblies and replicates - multiqc report\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/atac_seq/complex_samples.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featureCounts 16

  printf "\ncontrol and merging of tecnical replicates\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/atac_seq/control.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2 7

  printf "\ninput control different across same condition\n"
  seq2science run atac-seq -n --configfile tests/$WF/genrich_macs2.yaml --snakemakeOptions quiet=True config={samples:tests/atac_seq/complex_samples2.tsv,create_qc_report:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 genrich_pileup 4
  assert_rulecount $1 macs2_callpeak 4

  test_ran=1
fi

if [ $1 = "scatac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=scatac_seq

  printf "\nscatac-seq default\n"
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 create_SNAP_object 1

  printf "\nqc multiqc report\n"
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/samples.tsv,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 2  # None for sample and twice for trep

  printf "\nmultiple assemblies\n"
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/assemblies.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 bwa_mem2_index 2
  assert_rulecount $1 create_SNAP_object 2

  printf "\nmultiple assemblies - multiqc\n"
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/assemblies.tsv,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 4  # none for sample and twice for trep

  printf "\nmultiple replicates\n"
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/dag_sample.tsv,technical_replicates:merge} | tee tests/local_test_results/${1}_dag  # nothing to merge
  assert_rulecount $1 merge_replicates 0
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/replicates.tsv,technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 bwa_mem2 2
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/replicates.tsv,technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 2
  assert_rulecount $1 bwa_mem2 1

  printf "\nmultiple replicates - multiqc report\n"
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/replicates.tsv,technical_replicates:merge,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 2  # none for sample and twice for trep

  printf "\nmultiple assemblies and replicates\n"
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/dag_sample.tsv,technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 bwa_mem2 1
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/complex_samples.tsv,technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 bwa_mem2 4
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/complex_samples.tsv,technical_replicates:merge} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 3
  assert_rulecount $1 bwa_mem2 2

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  seq2science run scatac-seq -n --configfile tests/scatac_seq/default_config.yaml --snakemakeOptions quiet=True config={samples:tests/scatac_seq/complex_samples.tsv,technical_replicates:merge,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 3  # twice for PE trep, once for SE trep

  test_ran=1
fi

if [ $1 = "rna-seq" ]; then

  # RNA-seq workflow
  WF=rna_seq

  printf "\nrna-seq default\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,dexseq:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 gene_id2name 1
  assert_rulecount $1 htseq_count 1
  assert_rulecount $1 dexseq_count 1

  printf "\naligners\n"
#  seq2science run rna-seq -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,dexseq:True} | tee tests/local_test_results/${1}_dag  # default
#  assert_rulecount $1 star_align 1
#  assert_rulecount $1 dexseq_count 1
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:hisat2,dexseq:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 hisat2_align 1
  assert_rulecount $1 dexseq_count 1

  printf "\nquantifiers\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,quantifier:salmon,dexseq:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 1
  assert_rulecount $1 dexseq_count 1
  # seq2science run rna-seq -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,quantifier:htseq,dexseq:True} | tee tests/local_test_results/${1}_dag  # default
  # assert_rulecount $1 htseq_count 1
  # assert_rulecount $1 dexseq_count 1
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=False config={aligner:star,quantifier:featurecounts,dexseq:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 featurecounts 1
  assert_rulecount $1 dexseq_count 1

  printf "\ndecoy aware salmon index\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions quiet=True config={samples:tests/alignment/dag_sample.tsv,fastq_dir:tests/local_test_results/fastq} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 decoy_transcripts 1

  printf "\ntrackhub\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 0
  assert_rulecount $1 trackhub 1
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,quantifier:salmon,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 1
  assert_rulecount $1 trackhub 1

  printf "\nmultiqc report\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 4

  printf "\ncustom assembly\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,quantifier:salmon,custom_genome_extension:tests/local_test_results/ERCC92.fa,custom_annotation_extension:tests/local_test_results/ERCC92.gtf} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 extend_genome 1
  assert_rulecount $1 extend_genome_annotation 1
  assert_rulecount $1 salmon_index 1
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:star,quantifier:htseq,custom_genome_extension:tests/local_test_results/ERCC92.fa} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 extend_genome 1
  assert_rulecount $1 extend_genome_annotation 1
  assert_rulecount $1 star_index 1
  seq2science run rna-seq --skip-rerun -n --configfile tests/alignment/default_config.yaml --snakemakeOptions quiet=True config={aligner:hisat2,quantifier:htseq,custom_annotation_extension:tests/local_test_results/ERCC92.gtf} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 extend_genome 1
  assert_rulecount $1 extend_genome_annotation 1
  assert_rulecount $1 hisat2_splice_aware_index 1

  printf "\ndifferential expression analysis\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 htseq_count 10
  assert_rulecount $1 deseq2 1
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={quantifier:salmon,technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 10
  assert_rulecount $1 deseq2 1

  printf "\nmultiple assemblies - DEA\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:keep,samples:tests/rna_seq/complex_samples.tsv,dexseq:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 rename_sample 3
  assert_rulecount $1 star_index 2
  assert_rulecount $1 star_align 10
  assert_rulecount $1 dexseq_count 10
  assert_rulecount $1 deseq2 2  # TODO: it should only run once, since the contrast doesn't work in assembly2

  printf "\nmultiple assemblies - trackhubs\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:keep,samples:tests/rna_seq/complex_samples.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 trackhub 1

  printf "\nmultiple assemblies - multiqc\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:keep,samples:tests/rna_seq/complex_samples.tsv,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 24
  assert_rulecount $1 multiqc 2

  printf "\nmultiple replicates - DEA \n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:keep} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 htseq_count 10
  assert_rulecount $1 deseq2 1
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,dexseq:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 htseq_count 8
  assert_rulecount $1 dexseq_count 8
  assert_rulecount $1 deseq2 1

  printf "\nmultiple replicates - trackhubs\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 trackhub 1

  printf "\nmultiple replicates - multiqc\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 24
  assert_rulecount $1 htseq_count 8

  printf "\nmultiple assemblies and replicates - DEA \n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:keep,samples:tests/rna_seq/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 merge_replicates 0
  assert_rulecount $1 htseq_count 10
  assert_rulecount $1 deseq2 2  # TODO: it should only run once, since the contrast doesn't work in assembly2
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,samples:tests/rna_seq/complex_samples.tsv,dexseq:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 htseq_count 8
  assert_rulecount $1 dexseq_count 8
  assert_rulecount $1 count_matrix_DEXseq 2
  assert_rulecount $1 deseq2 2  # TODO: it should only run once, since the contrast doesn't work in assembly2
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={quantifier:salmon,technical_replicates:merge,samples:tests/rna_seq/complex_samples.tsv} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 8
  assert_rulecount $1 linked_txome 2
  assert_rulecount $1 deseq2 2  # TODO: it should only run once, since the contrast doesn't work in assembly2

  printf "\nmultiple assemblies and replicates - trackhub\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={quantifier:salmon,technical_replicates:merge,samples:tests/rna_seq/complex_samples.tsv,create_trackhub:True} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 salmon_quant 8
  assert_rulecount $1 star_align 8
  assert_rulecount $1 trackhub 1

  printf "\nmultiple assemblies and replicates - multiqc report\n"
  seq2science run rna-seq --skip-rerun -n --configfile tests/$WF/rna_seq_config.yaml --snakemakeOptions quiet=True config={technical_replicates:merge,samples:tests/rna_seq/complex_samples.tsv,create_qc_report:True,trimmer:trimgalore} | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastqc 24

  test_ran=1
fi

if [ $1 = "scrna-seq" ]; then

  # RNA-seq workflow
  WF=scrna_seq

  printf "\nscrna-seq default\n"
  seq2science run scrna-seq -n --configfile tests/scrna_seq/config.yaml --snakemakeOptions quiet=True | tee tests/local_test_results/${1}_dag
  assert_rulecount $1 fastp_SE 2
  assert_rulecount $1 fastq_pair 2
  assert_rulecount $1 kallistobus_ref 1
  assert_rulecount $1 kallistobus_count 2
  assert_rulecount $1 multiqc 1

  test_ran=1
fi

if [ $1 = "explain" ]; then
  seq2science init download-fastq --yes
  seq2science explain download-fastq

  seq2science init alignment --yes
  seq2science explain alignment

  seq2science init atac-seq --yes
  seq2science explain atac-seq

  seq2science init chip-seq --yes
  seq2science explain chip-seq

  seq2science init rna-seq --yes
  seq2science explain rna-seq

  seq2science init scrna-seq --yes
  seq2science explain scrna-seq

  seq2science init scatac-seq --yes
  seq2science explain scatac-seq

  test_ran=1
fi

# check if any test has run
if [ -z "$test_ran" ]; then
  printf "\nunrecognized input: ${1}\n"; exit
fi
