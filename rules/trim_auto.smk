# TODO: cutadapt/autodetect_adapter.py is a lot of code for what is does
rule trim_auto:
    input:
        download_fastq(expand("{result_dir}/fastq/{{sample}}.fastq.gz", **config))
    output:
        fastq=expand("{result_dir}/trimmed/{{sample}}_trimmed.fastq.gz", **config),
        qc=   expand("{result_dir}/trimmed/{{sample}}.fastq.gz_trimming_report.txt", **config)
    params:
        adapter=config['cut_adapter'],
        extra=  config['cut_params']
    log:
        expand("{log_dir}/trim_auto/{{sample}}.log", **config)
    threads: 6
    wrapper:
        "file:../../wrappers/cutadapt"