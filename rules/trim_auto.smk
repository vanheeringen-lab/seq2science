# TODO: cutadapt/autodetect_adapter.py is a lot of code for what is does

# checkpoint trim_auto_split:
#     input:
#         download_fastq(expand("{result_dir}/fastq/{{sample}}_{fqext1}_trimmed.{fqsuffix}.gz", **config),
#                        expand("{result_dir}/fastq/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config))
#     output:
#         fastq=expand("{result_dir}/trimmed/{{sample}}_trimmed.{fqsuffix}.gz", **config),
#         qc=   expand("{result_dir}/trimmed/{{sample}}.{fqsuffix}.gz_trimming_report.txt", **config)
#     params:
#         adapter=config['cut_adapter'],
#         extra=  config['cut_params']
#     log:
#         expand("{log_dir}/trim_auto/{{sample}}.log", **config)
#     threads: 6
#     wrapper:
#         "file:../../wrappers/cutadapt"


checkpoint trim_auto_splot:
    input:
        download_fastq(expand("{result_dir}/fastq/{{sample}}.{fqsuffix}.gz", **config))
    output:
        fastq=expand("{result_dir}/fastq/{{sample}}_trimmed.{fqsuffix}.gz", **config),
        qc=   expand("{result_dir}/trimmed/{{sample}}.{fqsuffix}.gz_trimming_report.txt", **config)
    params:
        adapter=config['cut_adapter'],
        extra=  config['cut_params']
    log:
        expand("{log_dir}/trim_auto/{{sample}}.log", **config)
    threads: 6
    wrapper:
        "file:../../wrappers/cutadapt"
