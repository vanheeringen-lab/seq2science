# to run these tests locally:
#   python setup.py develop
#   python setup.py build
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
  seq2science run download_fastq --cores $CORES --configfile tests/$WF/default_config.yaml

  WF=alignment

  # test genome & annotation downloading
  printf "\ndownload genome & annotation\n"
  seq2science run alignment --cores $CORES --configfile tests/$WF/remote_genome_n_sample.yaml --snakemakeOptions until=[get_genome]

  test_ran=1
fi

if [ $1 = "prep_align" ]; then

  WF=alignment

  printf "\ncreating environments\n"
  seq2science run $WF --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:bowtie2}
  seq2science run $WF --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:bwa-mem}
  seq2science run $WF --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:hisat2}
  seq2science run $WF --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:star}

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

  seq2science run $WF --cores 12 --configfile tests/$WF/default_config.yaml --snakemakeOptions overwrite_threads={bowtie2_align:$a,samtools_presort:$s} config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

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

  seq2science run $WF --cores 12 --configfile tests/$WF/default_config.yaml --snakemakeOptions overwrite_threads={bwa_mem:$a,samtools_presort:$s} config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

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

  seq2science run $WF --cores 12 --configfile tests/$WF/default_config.yaml --snakemakeOptions overwrite_threads={hisat2_align:$a,samtools_presort:$s} config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

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

  seq2science run $WF --cores 12 --configfile tests/$WF/default_config.yaml --snakemakeOptions overwrite_threads={star_align:$a,samtools_presort:$s} config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

  test_ran=1
fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq
  ALIGNER=bowtie2

  echo "Requires better test sample(s)"

#  printf "\natac-seq default\n"
#  skipped because the whole workflow needs to run again for multiqc
#  seq2science run atac_seq --cores $CORES --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2}

  printf "\natac-seq - multiqc report\n"
  seq2science run atac_seq --cores $CORES --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2,create_qc_report:True}

  printf "\natac-seq - trackhub\n"
  seq2science run atac_seq --cores $CORES --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2,create_trackhub:True}

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

# TODO: test samples are too similar for blind clustering and deseq2
  printf "\nrna-seq default - salmon\n"
  # TODO: currently omits blind clustering
  seq2science run rna_seq --cores $CORES --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={counts_dir:tests/local_test_results/salmon_counts} until=[txi_count_matrix]

#  printf "\nrna-seq default - salmon deseq2\n"
#  seq2science run rna_seq --cores $CORES --configfile tests/rna_seq/deseq2_config.yaml --snakemakeOptions config={counts_dir:tests/local_test_results/salmon_counts}

  # test STAR
  printf "\nrna-seq default - star\n"
  seq2science run rna_seq --cores $CORES --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={quantifier:star}

  printf "\nrna-seq default - star deseq2\n"
  seq2science run rna_seq --cores $CORES --configfile tests/rna_seq/deseq2_config.yaml --snakemakeOptions config={quantifier:star}

  printf "\nrna-seq default - trackhub\n"
  seq2science run rna_seq --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions config={quantifier:star,samples:tests/alignment/local_sample.tsv,genome_dir:tests,fastq_dir:../tinyfastq,create_trackhub:True}

  printf "\nrna-seq default - multiqc report\n"
  seq2science run rna_seq --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions config={quantifier:star,samples:tests/alignment/local_sample.tsv,genome_dir:tests,fastq_dir:../tinyfastq,create_qc_report:True}

  test_ran=1
fi

# check if any test has run
if [ -z "$test_ran" ]; then
  printf "\nunrecognized input: ${1}\n"; exit
fi
