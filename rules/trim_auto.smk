# TODO: cutadapt/autodetect_adapter.py is a lot of code for what is does
rule trim_auto:
    input:
        download_fastq(expand("{fastq_dir}/fastq/{{sample}}.fastq.gz", **config))
    output:
        fastq=expand("{fastq_dir}/trimmed/{{sample}}_trimmed.fastq.gz", **config),
        qc=   expand("{fastq_dir}/trimmed/{{sample}}.fastq.gz_trimming_report.txt", **config)
    params:
        adapter=config['cut_adapter'],
        extra=  config['cut_params']
    log:
        "logs/trim_auto/{sample}.log"
    threads: 6
    wrapper:
        "file:../../wrappers/cutadapt"