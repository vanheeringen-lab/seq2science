# to run these tests locally:
#   bash ./tests/dag_tests.sh TEST

# check if an argument was passed
if [ -z "$1" ]; then
    echo "No test specified"; exit
fi

CORES=28
set -e  # Exit immediately if a command exits with a non-zero status.

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
  test_ran=1
fi

if [ $1 = "cleanup_envs" ]; then
  rm -rf .snakemake
  rm -rf seq2science/workflows/download_fastq/.snakemake
  rm -rf seq2science/workflows/alignment/.snakemake
  rm -rf seq2science/workflows/atac_seq/.snakemake
  rm -rf seq2science/workflows/chip_seq/.snakemake
  rm -rf seq2science/workflows/rna_seq/.snakemake
  rm -rf seq2science/workflows/scATAC_seq/.snakemake
  test_ran=1
fi

if [ $1 = "download" ]; then

  WF=download_fastq

  # test basic downloading 1 PE and 1 SE
  printf "\ndownload SE and PE fastqs\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --configfile tests/$WF/default_config.yaml \
  --config samples=../../../tests/download_fastq/remote_samples.tsv \
  --use-conda --conda-frontend mamba -j $CORES

  WF=alignment

  # test genome & annotation downloading
  printf "\ndownload genome & annotation\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --configfile tests/$WF/remote_genome_n_sample.yaml \
  --use-conda --conda-frontend mamba -j $CORES \
  --until get_genome

  test_ran=1
fi

if [ $1 = "prep_align" ]; then

  WF=alignment

  printf "\ncreating environments\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  -j $CORES --quiet --configfile tests/$WF/default_config.yaml \
  --config samples=../../../tests/$WF/remote_genome_n_sample.tsv aligner=bowtie2 \
  --use-conda --conda-create-envs-only --conda-frontend mamba

  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  -j $CORES --quiet --configfile tests/$WF/default_config.yaml \
  --config samples=../../../tests/$WF/remote_genome_n_sample.tsv aligner=bwa-mem \
  --use-conda --conda-create-envs-only --conda-frontend mamba

  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  -j $CORES --quiet --configfile tests/$WF/default_config.yaml \
  --config samples=../../../tests/$WF/remote_genome_n_sample.tsv aligner=hisat2 \
  --use-conda --conda-create-envs-only --conda-frontend mamba

  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  -j $CORES --quiet --configfile tests/$WF/default_config.yaml \
  --config samples=../../../tests/$WF/remote_genome_n_sample.tsv aligner=star \
  --use-conda --conda-create-envs-only --conda-frontend mamba

  test_ran=1
fi

if [ $1 = "bowtie2" ]; then

  ALIGNER=bowtie2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
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
      fastq_dir=../../../tests/tinyfastq \
      genome_dir=../../../tests \
      result_dir=../../../${RESULTS_DIR} \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

  test_ran=1
fi

if [ $1 = "bwa-mem" ]; then

  ALIGNER=bwa-mem
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
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
#      fastq_dir=../../../tests/tinyfastq \
#      genome_dir=../../../tests \
#      result_dir=../../../${RESULTS_DIR} \
#  -j $c --set-threads bwa_mem=$a samtools_presort=$s

  test_ran=1
fi

if [ $1 = "hisat2" ]; then

  ALIGNER=hisat2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
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
      fastq_dir=../../../tests/tinyfastq \
      genome_dir=../../../tests \
      result_dir=../../../${RESULTS_DIR} \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

  test_ran=1
fi

if [ $1 = "star" ]; then

  ALIGNER=star
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
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
      fastq_dir=../../../tests/tinyfastq \
      genome_dir=../../../tests \
      result_dir=../../../${RESULTS_DIR} \
  -j $c --set-threads ${ALIGNER}_align=$a samtools_presort=$s

  test_ran=1
fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq
  ALIGNER=bowtie2

  echo "Requires better test sample(s)"

#  skipped because the whole workflow needs to run again for multiqc
#  printf "\natac-seq default\n"
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda --conda-frontend mamba -j $CORES \
#  --configfile \
#      tests/alignment/remote_genome_n_sample.yaml \
#  --config \
#      aligner=bowtie2

# TODO: error!
#  printf "\natac-seq - multiqc report\n"
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda --conda-frontend mamba -j $CORES \
#  --configfile \
#      tests/alignment/remote_genome_n_sample.yaml \
#  --config \
#      aligner=bowtie2 \
#      create_qc_report=True

# TODO: error!
#  printf "\natac-seq - trackhub\n"
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda --conda-frontend mamba -j $CORES \
#  --configfile \
#      tests/alignment/remote_genome_n_sample.yaml \
#  --config \
#      aligner=bowtie2 \
#      create_trackhub=True

  test_ran=1
fi

if [ $1 = "scatac-seq" ]; then

  # scATAC-seq workflow
  WF=scATAC_seq

  echo "Requires test sample(s)"

  test_ran=1
fi

if [ $1 = "rna-seq" ]; then

  # RNA-seq workflow
  WF=rna_seq

  printf "\nrna-seq default - salmon\n"
# test samples are too similar for deseq2
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda --conda-frontend mamba -j $CORES \
#  --configfile \
#      tests/rna_seq/salmon_config.yaml \
#  --config \
#      counts_dir=$(pwd)/tests/local_test_results/salmon_counts

# TODO: error: ChildIOException: tinydata.gtf linked_txome/get_annotation
#  printf "\nrna-seq default - salmon deseq2\n"
#  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
#  --use-conda --conda-frontend mamba -j $CORES \
#  --configfile \
#      tests/rna_seq/deseq2_config.yaml \
#  --config \
#      counts_dir=$(pwd)/tests/local_test_results/salmon_counts

  # test STAR
  printf "\nrna-seq default - star\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda --conda-frontend mamba -j $CORES \
  --configfile \
      tests/rna_seq/salmon_config.yaml \
  --config quantifier=star

  printf "\nrna-seq default - star deseq2\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda --conda-frontend mamba -j $CORES \
  --configfile \
      tests/rna_seq/deseq2_config.yaml \
  --config quantifier=star

  printf "\nrna-seq default - trackhub\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda --conda-frontend mamba -j $CORES \
  --configfile \
      tests/alignment/default_config.yaml \
  --config \
    samples=../../../tests/alignment/local_sample.tsv \
    genome_dir=../../../tests \
    fastq_dir=../../tests/tinyfastq \
    create_trackhub=True

  printf "\nrna-seq default - multiqc report\n"
  snakemake -s seq2science/workflows/$WF/Snakefile --directory seq2science/workflows/$WF \
  --use-conda --conda-frontend mamba -j $CORES \
  --configfile \
      tests/alignment/default_config.yaml \
  --config \
    samples=../../../tests/alignment/local_sample.tsv \
    genome_dir=../../../tests \
    fastq_dir=../../tests/tinyfastq \
    create_qc_report=True

  test_ran=1
fi

# check if any test has run
if [ -z "$test_ran" ]; then
  printf "\nunrecognized input: ${1}\n"; exit
fi
