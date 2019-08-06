rule samtools_stats:
    input:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.bam", **config)
    output:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.samtools_stats.txt", **config)
    log:
        expand("{log_dir}/samtools_stats/{{sample}}-{{assembly}}.log", **config)
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"


rule featureCounts:
    # https://www.biostars.org/p/337872/
    # https://www.biostars.org/p/228636/
    input:
        bam=expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.bam", **config),
        peak=expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}_peaks.narrowPeak", **config)
    output:
        tmp_saf=temp(expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}.saf", **config)),
        real_out=expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}_featureCounts.txt", **config),
        summary=expand("{result_dir}/{{peak_caller}}/{{sample}}-{{assembly}}_featureCounts.txt.summary", **config)
    log:
        expand("{log_dir}/featureCounts/{{sample}}-{{assembly}}-{{peak_caller}}.log", **config)
    threads: 4
    conda:
        "../envs/subread.yaml"
    shell:
        """
        ## Make a custom "SAF" file which featureCounts needs:
        awk 'BEGIN{{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' {input.peak} 1> {output.tmp_saf} 2> {log}

        ## run featureCounts
        featureCounts -T {threads} -p -a {output.tmp_saf} -F SAF -o {output.real_out} {input.bam} > {log} 2>&1
        """


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
