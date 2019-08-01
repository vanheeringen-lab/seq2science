rule fastqc:
    input:
        f"{{path}}/{{fname}}.{config['fqsuffix']}.gz"
    output:
        "{path}/{fname}_fastqc.html",
        "{path}/{fname}_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    shell:
        "fastqc {input} -O {wildcards.path} --quiet"


def get_qc_files(wildcards):
    assert 'quality_control' in globals(), "When trying to generate multiqc output, make sure that the "\
                                           "variable 'quality_control' exists and contains all the "\
                                           "relevant quality control files."
    return quality_control


rule multiqc:
    input:
       get_qc_files
    output:
        expand("{result_dir}/qc/multiqc.html", **config),
        directory(expand("{result_dir}/qc/multiqc_data", **config))
    params:
        "{result_dir}/qc/".format(**config)
    log:
        expand("{log_dir}/multiqc.log", **config)
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc {input} -o {params} -n multiqc.html --config ../../multiqc_config.yaml > {log} 2>&1"


def get_trimming_qc(sample):
    if config['layout'][sample] == 'SINGLE':
        return expand([f"{{result_dir}}/{{trimmed_dir}}/{sample}_fastqc.zip",
                       f"{{result_dir}}/{{trimmed_dir}}/{sample}_trimmed_fastqc.zip",
                       f"{{result_dir}}/{{trimmed_dir}}/{sample}.{{fqsuffix}}.gz_trimming_report.txt"],
                       **config)
    else:
        return expand([f"{{result_dir}}/{{trimmed_dir}}/{sample}_{{fqext}}_fastqc.zip",
                       f"{{result_dir}}/{{trimmed_dir}}/{sample}_{{fqext}}_trimmed_fastqc.zip",
                       f"{{result_dir}}/{{trimmed_dir}}/{sample}_{{fqext}}.{{fqsuffix}}.gz_trimming_report.txt"],
                       **config)

def get_alignment_qc(sample):
    return expand([f"{{result_dir}}/{{dedup_dir}}/{sample}-{samples.loc[sample]['assembly']}.metrics.txt",
                   f"{{result_dir}}/{{dedup_dir}}/{sample}-{samples.loc[sample]['assembly']}.samtools_stats.txt"],
                   **config)


def get_peak_calling_qc(sample):
    return expand(f"{{result_dir}}/{{peak_caller}}/{sample}-{samples.loc[sample]['assembly']}_featureCounts.txt.summary",
                  **config)
