# #!/usr/bin/env bash

# to run these tests locally:
# 1)   chmod u+x ./tests/dag_tests.sh
# 2)   ./tests/dag_tests.sh TEST

if [ -z "$1" ]
  then
    echo "No test specified"
    exit
fi

CORES=48
set -e  # Exit immediately if a command exits with a non-zero status.
#function assert_rulecount {
#  # check if the DAG (stored with  | tee tests/local_test_results/val  ) ran rule $1 exactly $2 times
#  val=$(cat tests/local_test_results/val | grep -w $1 | cut -f2);
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
#mkdir -p tests/local_test_results

if [ $1 = "cleanup_files" ]; then
  rm -rf tests/local_test_results
  rm -rf tests/tinydata/index
  rm -rf tests/tinydata/decoy_transcripts
  rm -rf tests/tinydata/tinydata.2bit
  rm -rf tests/tinydata/cytoBandIdeo.bb
  rm -rf tests/tinydata/tinydata.bb
  rm -rf tests/tinydata/tinydata.custom*
  rm -rf tests/tinydata/tinydata.gc5Base.bw
  rm -rf tests/tinydata/tinydata.gtf
  rm -rf tests/tinydata/tinydata.ix
  rm -rf tests/tinydata/tinydata.ixx
  rm -rf tests/tinydata/tinydata.transcripts.fa
  rm -rf tests/tinydata/tinydata_softmasking.bb
fi

if [ $1 = "cleanup_envs" ]; then
  rm -rf .snakemake
  rm -rf seq2science/workflows/download_fastq/.snakemake
  rm -rf seq2science/workflows/alignment/.snakemake
  rm -rf seq2science/workflows/atac_seq/.snakemake
  rm -rf seq2science/workflows/chip_seq/.snakemake
  rm -rf seq2science/workflows/rna_seq/.snakemake
  rm -rf seq2science/workflows/scATAC_seq/.snakemake
fi

if [ $1 = "download" ]; then

  WF=download_fastq

  # test basic downloading 1 PE and 1 SE
  printf "\ndownload SE and PE fastqs\n\n"
  snakemake --use-conda -j $CORES -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --configfile tests/$WF/default_config.yaml \
  --config samples=../../../tests/download_fastq/remote_samples.tsv

  WF=alignment

  # test genome & annotation downloading
  printf "\ndownload genome & annotation\n\n"
  snakemake --use-conda -j $CORES -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --configfile tests/$WF/remote_genome_n_sample.yaml \
  --until get_genome

fi

if [ $1 = "prep_align" ]; then

  snakemake -s seq2science/workflows/alignment/Snakefile --directory seq2science/workflows/alignment \
  --use-conda -j $CORES --configfile tests/alignment/default_config.yaml \
  --config samples=../../../tests/alignment/remote_genome_n_sample.tsv aligner=bowtie2 \
  --conda-create-envs-only --quiet

  snakemake -s seq2science/workflows/alignment/Snakefile --directory seq2science/workflows/alignment \
  --use-conda -j $CORES --configfile tests/alignment/default_config.yaml \
  --config samples=../../../tests/alignment/remote_genome_n_sample.tsv aligner=bwa-mem \
  --conda-create-envs-only --quiet

  snakemake -s seq2science/workflows/alignment/Snakefile --directory seq2science/workflows/alignment \
  --use-conda -j $CORES --configfile tests/alignment/default_config.yaml \
  --config samples=../../../tests/alignment/remote_genome_n_sample.tsv aligner=hisat2 \
  --conda-create-envs-only --quiet

  snakemake -s seq2science/workflows/alignment/Snakefile --directory seq2science/workflows/alignment \
  --use-conda -j $CORES --configfile tests/alignment/default_config.yaml \
  --config samples=../../../tests/alignment/remote_genome_n_sample.tsv aligner=star \
  --conda-create-envs-only --quiet

fi

if [ $1 = "bowtie2" ]; then

  ALIGNER=bowtie2
  WF=alignment
  let "c = $CORES / 4"
  let "a = $c - 2"
  let "s = 2"

  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda --nolock --notemp \
  --configfile \
      tests/$WF/default_config.yaml \
  --config \
      aligner=$ALIGNER \
      samples=../../../tests/alignment/local_sample.tsv \
      fastq_dir=$(pwd)/tests/tinyfastq \
      genome_dir=$(pwd)/tests \
      result_dir=$(pwd)/tests/local_test_results/$ALIGNER \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

fi

if [ $1 = "bwa-mem" ]; then

  ALIGNER=bwa
  WF=alignment
  let "c = $CORES / 4"
  let "a = $c - 2"
  let "s = 2"

  # TODO: error!
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda --nolock --notemp \
#  --configfile \
#      tests/$WF/default_config.yaml \
#  --config \
#      aligner=$ALIGNER \
#      samples=../../../tests/alignment/local_sample.tsv \
#      fastq_dir=$(pwd)/tests/tinyfastq \
#      genome_dir=$(pwd)/tests \
#      result_dir=$(pwd)/tests/local_test_results/$ALIGNER \
#  -j $c --set-threads bwa_mem=$a samtools_presort=$s

fi

if [ $1 = "hisat2" ]; then

  ALIGNER=hisat2
  WF=alignment
  let "c = $CORES / 4"
  let "a = $c - 2"
  let "s = 2"

  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda --nolock --notemp \
  --configfile \
      tests/$WF/default_config.yaml \
  --config \
      aligner=$ALIGNER \
      samples=../../../tests/alignment/local_sample.tsv \
      fastq_dir=$(pwd)/tests/tinyfastq \
      genome_dir=$(pwd)/tests \
      result_dir=$(pwd)/tests/local_test_results/$ALIGNER \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

fi

if [ $1 = "star" ]; then

  ALIGNER=star
  WF=alignment
  let "c = $CORES / 4"
  let "a = $c - 2"
  let "s = 2"

  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda --nolock --notemp \
  --configfile \
      tests/$WF/default_config.yaml \
  --config \
      aligner=$ALIGNER \
      samples=../../../tests/alignment/local_sample.tsv \
      fastq_dir=$(pwd)/tests/tinyfastq \
      genome_dir=$(pwd)/tests \
      result_dir=$(pwd)/tests/local_test_results/$ALIGNER \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq
  ALIGNER=bowtie2

  echo "Requires better test sample(s)"

#  skipped because the whole workflow needs to run again for multiqc
#  printf "\natac-seq default\n"
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda -j $CORES \
#  --configfile \
#      tests/alignment/remote_genome_n_sample.yaml \
#  --config \
#      aligner=bowtie2

  printf "\natac-seq - multiqc report\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda -j $CORES \
  --configfile \
      tests/alignment/remote_genome_n_sample.yaml \
  --config \
      aligner=bowtie2 \
      create_qc_report=True

  printf "\natac-seq - trackhub\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda -j $CORES \
  --configfile \
      tests/alignment/remote_genome_n_sample.yaml \
  --config \
      aligner=bowtie2 \
      create_trackhub=True

fi

if [ $1 = "scatac-seq" ]; then

  # scATAC-seq workflow
  WF=scATAC_seq

  echo "Requires test sample(s)"

fi

if [ $1 = "rna-seq" ]; then

  # RNA-seq workflow
  WF=rna_seq

  printf "\nrna-seq default - salmon\n"
# TODO: error!
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda -j $CORES \
#  --configfile \
#      tests/rna_seq/salmon_config.yaml \
#  --omit-from blind_clustering deseq2

#  test samples are too similar for deseq2
#  printf "\nrna-seq default - salmon deseq2\n"
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda -j $CORES \
#  --configfile \
#      tests/rna_seq/salmon_config.yaml

  rm -rf tests/local_test_results/gene_counts

  # test STAR
  printf "\nrna-seq default - star\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda -j $CORES \
  --configfile \
      tests/rna_seq/salmon_config.yaml \
  --config quantifier=star \
  --omit-from blind_clustering deseq2

  printf "\nrna-seq default - star deseq2\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda -j $CORES \
  --configfile \
      tests/rna_seq/salmon_config.yaml \
  --config quantifier=star

  printf "\nrna-seq default - trackhub\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda -j $CORES \
  --configfile \
      tests/alignment/default_config.yaml \
  --config \
    samples=../../../tests/alignment/local_sample.tsv \
    genome_dir=../../tests \
    fastq_dir=../tests/tinyfastq \
    create_trackhub=True

  printf "\nrna-seq default - multiqc report\n"
  # TODO: error!
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda -j $CORES \
#  --configfile \
#      tests/alignment/default_config.yaml \
#  --config \
#    samples=../../../tests/alignment/local_sample.tsv \
#    genome_dir=../../tests \
#    fastq_dir=../tests/tinyfastq \
#    create_qc_report=True

fi
