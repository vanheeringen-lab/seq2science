# to run these tests locally:
#   python setup.py develop
#   bash ./tests/run_tests.sh TEST

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
  rm -rf tests/tinydata/cytoBandIdeo.bb
  rm -rf tests/tinydata/tinydata.2bit
  rm -rf tests/tinydata/tinydata.bb
  rm -rf tests/tinydata/tinydata.custom*
  rm -rf tests/tinydata/tinydata.gc5Base.bw
  rm -rf tests/tinydata/tinydata.ix
  rm -rf tests/tinydata/tinydata.ixx
  rm -rf tests/tinydata/tinydata.transcripts.fa
  rm -rf tests/tinydata/tinydata_softmasking.bb
  test_ran=1
fi

if [ $1 = "cleanup_envs" ]; then
  rm -rf .snakemake
  rm -rf seq2science/.snakemake
  test_ran=1
fi

if [ $1 = "download" ]; then

  WF=download_fastq

  # test basic downloading 1 PE and 1 SE
  printf "\ndownload SE and PE fastqs\n"
  seq2science run download-fastq --cores $CORES --configfile tests/$WF/default_config.yaml

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
  seq2science run $WF --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:bwa-mem2}
  seq2science run $WF --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:hisat2}
  seq2science run $WF --cores $CORES --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:star}

  test_ran=1
fi

if [ $1 = "bowtie2" ]; then

  ALIGNER=bowtie2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 5"

  seq2science run $WF --cores $c --configfile tests/$WF/default_config.yaml --snakemakeOptions config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

  test_ran=1
fi

if [ $1 = "bwa-mem1" ]; then

  ALIGNER=bwa-mem
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 5"

  seq2science run $WF --cores $c --configfile tests/$WF/default_config.yaml --snakemakeOptions config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

  test_ran=1
fi

if [ $1 = "bwa-mem2" ]; then

  ALIGNER=bwa-mem2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 5"

  seq2science run $WF --cores 12 --configfile tests/$WF/default_config.yaml --snakemakeOptions resources={mem_gb:999} config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

  test_ran=1
fi

if [ $1 = "hisat2" ]; then

  ALIGNER=hisat2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 5"

  seq2science run $WF --cores $c --configfile tests/$WF/default_config.yaml --snakemakeOptions config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

  test_ran=1
fi

if [ $1 = "star" ]; then

  ALIGNER=star
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 5"

  seq2science run $WF --cores $c --configfile tests/$WF/default_config.yaml --snakemakeOptions config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../../tests/tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER}

  test_ran=1
fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq
  ALIGNER=bowtie2

  echo "Requires better test sample(s)"

#  printf "\natac-seq default\n"
#  skipped because the whole workflow needs to run again for multiqc
#  seq2science run atac-seq --cores $CORES --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2}

  printf "\natac-seq - multiqc report\n"
  seq2science run atac-seq --cores $CORES --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2,create_qc_report:True}

  printf "\natac-seq - trackhub\n"
  seq2science run atac-seq --cores $CORES --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2,create_trackhub:True}

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
  printf "\nrna-seq default - quantification\n"
  seq2science run rna-seq --cores $CORES --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={counts_dir:salmon_counts} \
  omit_from=[blind_clustering]  # <- remove when fixed

  printf "\nrna-seq default - quantification deseq2\n"
  seq2science run rna-seq --cores $CORES --configfile tests/rna_seq/deseq2_config.yaml --snakemakeOptions config={counts_dir:salmon_counts} \
  omit_from=[blind_clustering,deseq2]  # <- remove when fixed

  printf "\nrna-seq default - counting\n"
  seq2science run rna-seq --cores $CORES --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={quantifier:htseq} \
  omit_from=[blind_clustering]  # <- remove when fixed
  seq2science run rna-seq --cores $CORES --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={quantifier:featurecounts,counts_dir:fc_counts} \
  omit_from=[blind_clustering]  # <- remove when fixed
  seq2science run rna-seq --cores $CORES --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={aligner:hisat2,quantifier:htseq,final_bam_dir:hisat2_final_bam,counts_dir:hisat2_counts} \
  omit_from=[blind_clustering]  # <- remove when fixed

  printf "\nrna-seq default - counting deseq2\n"
  seq2science run rna-seq --cores $CORES --configfile tests/rna_seq/deseq2_config.yaml \
  --snakemakeOptions omit_from=[blind_clustering,deseq2]  # <- remove when fixed

  printf "\nrna-seq default - trackhub\n"
  seq2science run rna-seq --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions config={samples:tests/alignment/stranded_sample.tsv,genome_dir:tests,fastq_dir:../tinyfastq,aligner:star,create_trackhub:True}

  printf "\nrna-seq default - multiqc report\n"
  seq2science run rna-seq --cores $CORES --configfile tests/alignment/default_config.yaml --snakemakeOptions config={samples:tests/alignment/local_sample.tsv,genome_dir:tests,fastq_dir:../tinyfastq,aligner:star,create_qc_report:True}

  test_ran=1
fi

# check if any test has run
if [ -z "$test_ran" ]; then
  printf "\nunrecognized input: ${1}\n"; exit
fi
