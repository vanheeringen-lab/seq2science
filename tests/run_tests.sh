# to run these tests locally:
#   mamba install toml --yes
#   python setup.py develop
#   bash ./tests/run_tests.sh TEST

# for lazy asses:
# seq2science clean && bash ./tests/run_tests.sh download && bash ./tests/run_tests.sh prep_align && bash ./tests/run_tests.sh bowtie2 && bash ./tests/run_tests.sh bwa-mem1 && bash ./tests/run_tests.sh bwa-mem2 && bash ./tests/run_tests.sh hisat2 && bash ./tests/run_tests.sh star && bash ./tests/run_tests.sh atac-seq && bash ./tests/run_tests.sh scatac-seq && bash ./tests/run_tests.sh rna-seq

# check if an argument was passed
if [ -z "$1" ]; then
    echo "No test specified"; exit
fi

CORES=28
set -e  # Exit immediately if a command exits with a non-zero status.
touch -m tests/tinydata/*  # update testdata timestamps (prevent unexpected reruns)

if [ $1 = "cleanup_files" ]; then
  rm -rf tests/local_test_results
  rm -rf tests/tinydata/decoy_transcripts
  rm -rf tests/tinydata/index
  rm -rf tests/tinydata/annotation.bigBed
  rm -rf tests/tinydata/annotation.ix
  rm -rf tests/tinydata/annotation.ixx
  rm -rf tests/tinydata/cytoBandIdeo.bb
  rm -rf tests/tinydata/gene_id2name.tsv
  rm -rf tests/tinydata/tinydata.2bit
  rm -rf tests/tinydata/tinydata.DEXseq*
  rm -rf tests/tinydata/tinydata.gc5Base.bw
  rm -rf tests/tinydata/tinydata.transcripts.fa
  rm -rf tests/tinydata/tinydata_softmasking.bb
  rm -rf tests/tinydata_custom
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
  seq2science run download-fastq --cores $CORES -r --configfile tests/$WF/default_config.yaml --snakemakeOptions show_failed_logs=True printshellcmds=True

  WF=alignment

  # test genome & annotation downloading
  printf "\ndownload genome & annotation\n"
  seq2science run alignment --cores $CORES -r --configfile tests/$WF/remote_genome_n_sample.yaml --snakemakeOptions until=[get_genome_annotation] show_failed_logs=True

  test_ran=1
fi

if [ $1 = "prep_align" ]; then

  WF=alignment

  printf "\ncreating environments\n"
  seq2science run $WF --cores $CORES -r --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:bowtie2}
  seq2science run $WF --cores $CORES -r --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:bwa-mem}
  seq2science run $WF --cores $CORES -r --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:bwa-mem2}
  seq2science run $WF --cores $CORES -r --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:chromap}
  seq2science run $WF --cores $CORES -r --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:hisat2}
  seq2science run $WF --cores $CORES -r --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:minimap2}
  seq2science run $WF --cores $CORES -r --configfile tests/$WF/default_config.yaml --snakemakeOptions quiet=True conda_create_envs_only=true config={samples:tests/$WF/remote_genome_n_sample.tsv,aligner:star}

  test_ran=1
fi

if [ $1 = "bowtie2" ]; then

  ALIGNER=bowtie2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 7"

  seq2science run $WF --cores $c -r --configfile tests/$WF/default_config.yaml --snakemakeOptions config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER} show_failed_logs=True

  test_ran=1
fi

if [ $1 = "bwa-mem1" ]; then

  ALIGNER=bwa-mem
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 7"

  seq2science run $WF --cores $c -r --configfile tests/$WF/default_config.yaml --snakemakeOptions until=[samtools_presort] config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER} show_failed_logs=True

  test_ran=1
fi

if [ $1 = "bwa-mem2" ]; then

  ALIGNER=bwa-mem2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 7"

  seq2science run $WF --cores $c -r --configfile tests/$WF/default_config.yaml --snakemakeOptions until=[samtools_presort] resources={mem_gb:999} config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER} show_failed_logs=True

  test_ran=1
fi

if [ $1 = "bwa-mem2" ]; then

  ALIGNER=chromap
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 7"

  seq2science run $WF --cores $c -r --configfile tests/$WF/default_config.yaml --snakemakeOptions until=[samtools_presort] resources={mem_gb:999} config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER} show_failed_logs=True

  test_ran=1
fi

if [ $1 = "hisat2" ]; then

  ALIGNER=hisat2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 7"

  seq2science run $WF --cores $c -r --configfile tests/$WF/default_config.yaml --snakemakeOptions until=[samtools_presort] config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER} show_failed_logs=True

  test_ran=1
fi

if [ $1 = "minimap2" ]; then

  ALIGNER=minimap2
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 7"

  seq2science run $WF --cores $c -r --configfile tests/$WF/default_config.yaml --snakemakeOptions until=[samtools_presort] resources={mem_gb:999} config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER} show_failed_logs=True

  test_ran=1
fi

if [ $1 = "star" ]; then

  ALIGNER=star
  WF=alignment
  RESULTS_DIR=tests/local_test_results/${ALIGNER}
  mkdir -p $RESULTS_DIR
  let "c = $CORES / 7"

  seq2science run $WF --cores $c -r --configfile tests/$WF/default_config.yaml --snakemakeOptions until=[samtools_presort] config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../tinyfastq,genome_dir:tests,result_dir:$RESULTS_DIR,aligner:$ALIGNER} show_failed_logs=True
  # seq2science run alignment -nr --configfile tests/alignment/default_config.yaml --snakemakeOptions until=[samtools_presort] config={samples:tests/alignment/local_sample.tsv,fastq_dir:../../tinyfastq,genome_dir:tests,result_dir:tests/local_test_results/star,aligner:star}

  test_ran=1
fi

if [ $1 = "atac-seq" ]; then

  # ATAC-seq workflow (also covers ChIP-seq workflow)
  WF=atac_seq
  ALIGNER=bowtie2

  echo "Requires better test sample(s)"

#  printf "\natac-seq default\n"
#  skipped because the whole workflow needs to run again for multiqc
#  seq2science run atac-seq --cores $CORES -r --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2}

  printf "\natac-seq - multiqc report\n"
  seq2science run atac-seq --cores $CORES -r --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2,create_qc_report:True,custom_genome_extension:tests/tinydata/tinyERCC92.fa,min_template_length:50} show_failed_logs=True

  printf "\natac-seq - trackhub\n"
  seq2science run atac-seq --cores $CORES -r --configfile tests/alignment/remote_genome_n_sample.yaml --snakemakeOptions config={aligner:bowtie2,create_trackhub:True,custom_genome_extension:tests/tinydata/tinyERCC92.fa,max_template_length:150} show_failed_logs=True

  # fake input data & fix timestamps
  mkdir -p tests/local_test_results/GRCh38.p13
  touch tests/local_test_results/GRCh38.p13/GRCh38.p13.fa
  touch tests/local_test_results/GRCh38.p13/GRCh38.p13.annotation.gtf
  touch tests/local_test_results/GRCh38.p13/GRCh38.p13.annotation.bed
  mkdir -p tests/local_test_results/fastq
  for ext in R1 R2; do
    for n in 1 2 3 4 5 6 7 8; do
      touch tests/local_test_results/fastq/sample${n}_${ext}.fastq.gz
    done
    touch tests/local_test_results/fastq/replicate_${ext}.fastq.gz
  done
  touch -m tests/deseq2/atac/macs2/GRCh38.p13_raw.tsv  # update timestamps
  touch -m tests/deseq2/atac/macs2/GRCh38.p13_onehotpeaks.tsv
  touch -m tests/deseq2/atac/macs2/GRCh38.p13_meancenter_log2_quantilenorm.tsv
  touch -m tests/deseq2/atac/macs2/GRCh38.p13_meancenter_log2_TMM.tsv
  touch -m tests/deseq2/atac/macs2/GRCh38.p13_meancenter_log2_RLE.tsv
  touch -m tests/deseq2/atac/macs2/GRCh38.p13_meancenter_log2_upperquartile.tsv

  # run DESeq2
  seq2science run atac-seq --skip-rerun --cores $CORES -r --configfile tests/deseq2/atac/config.yaml --snakemakeOptions config={deseq2_dir:deseq_atac} show_failed_logs=True until=[deseq2]

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
  seq2science run rna-seq --cores $CORES -r --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={counts_dir:salmon_counts,custom_genome_extension:tests/tinydata/tinyERCC92.fa,custom_annotation_extension:tests/tinydata/tinyERCC92.gtf} show_failed_logs=True \
  omit_from=[blind_clustering]  # <- remove when fixed

  printf "\nrna-seq default - quantification deseq2\n"
#  seq2science run rna-seq --cores $CORES -r --configfile tests/rna_seq/deseq2_config.yaml --snakemakeOptions config={counts_dir:salmon_counts} show_failed_logs=True

  # fake input data & fix timestamps
  mkdir -p tests/local_test_results/GRCh38.p13
  touch tests/local_test_results/GRCh38.p13/GRCh38.p13.fa
  touch tests/local_test_results/GRCh38.p13/GRCh38.p13.annotation.gtf
  touch tests/local_test_results/GRCh38.p13/GRCh38.p13.annotation.bed
  touch tests/local_test_results/GRCh38.p13/gene_id2name.tsv
  mkdir -p tests/local_test_results/fastq
  for ext in R1 R2; do
    for n in 1 2 3 4 5 6 7 8; do
      touch tests/local_test_results/fastq/sample${n}_${ext}.fastq.gz
    done
    touch tests/local_test_results/fastq/replicate_${ext}.fastq.gz
  done
  touch -m tests/deseq2/rna/counts/GRCh38.p13-counts.tsv  # update timestamp

  # run blind clustering
  seq2science run rna-seq --skip-rerun --cores $CORES -r --configfile tests/deseq2/rna/config.yaml --snakemakeOptions config={deseq2_dir:deseq_rna,create_qc_report:true} show_failed_logs=True \
  targets=[$(pwd)/tests/local_test_results/qc/plotCorrelation/GRCh38.p13-DESeq2_sample_distance_clustering_mqc.png]

  # run DESeq2
  seq2science run rna-seq --skip-rerun --cores $CORES -r --configfile tests/deseq2/rna/config.yaml --snakemakeOptions config={deseq2_dir:deseq_rna} show_failed_logs=True

  printf "\nrna-seq default - counting\n"
  seq2science run rna-seq --cores $CORES -r --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={quantifier:htseq,dexseq:True,custom_genome_extension:tests/tinydata/tinyERCC92.fa,custom_annotation_extension:tests/tinydata/tinyERCC92.gtf} show_failed_logs=True \
  omit_from=[blind_clustering]  # <- remove when fixed
  seq2science run rna-seq --cores $CORES -r --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={quantifier:featurecounts,counts_dir:fc_counts} show_failed_logs=True \
  omit_from=[blind_clustering]  # <- remove when fixed
  seq2science run rna-seq --cores $CORES -r --configfile tests/rna_seq/salmon_config.yaml --snakemakeOptions config={aligner:hisat2,quantifier:htseq,final_bam_dir:hisat2_final_bam,counts_dir:hisat2_counts,custom_genome_extension:tests/tinydata/tinyERCC92.fa,custom_annotation_extension:tests/tinydata/tinyERCC92.gtf} show_failed_logs=True \
  omit_from=[blind_clustering]  # <- remove when fixed

  printf "\nrna-seq default - counting deseq2\n"
  seq2science run rna-seq --cores $CORES -r --configfile tests/rna_seq/deseq2_config.yaml --snakemakeOptions show_failed_logs=True \
  omit_from=[blind_clustering,deseq2]  # <- remove when fixed

  printf "\nrna-seq default - trackhub\n"
  # deeptools: removed normalization due to test samples being too small
  seq2science run rna-seq --cores $CORES -r --configfile tests/alignment/default_config.yaml --snakemakeOptions config={samples:tests/alignment/stranded_sample.tsv,genome_dir:tests,fastq_dir:../tinyfastq,aligner:star,create_trackhub:True,custom_genome_extension:tests/tinydata/tinyERCC92.fa,custom_annotation_extension:tests/tinydata/tinyERCC92.gtf,deeptools_flags:-v} show_failed_logs=True

  printf "\nrna-seq default - multiqc report\n"
  seq2science run rna-seq --cores $CORES -r --configfile tests/alignment/default_config.yaml --snakemakeOptions config={samples:tests/alignment/local_sample.tsv,genome_dir:tests,fastq_dir:../tinyfastq,aligner:star,create_qc_report:True,custom_genome_extension:tests/tinydata/tinyERCC92.fa,custom_annotation_extension:tests/tinydata/tinyERCC92.gtf} show_failed_logs=True

  test_ran=1
fi

# check if any test has run
if [ -z "$test_ran" ]; then
  printf "\nunrecognized input: ${1}\n"; exit
fi
