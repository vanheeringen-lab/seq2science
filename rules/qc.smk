rule samtools_stats:
    """
    Get general stats from bam files like percentage mapped.
    """
    input:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{{bam_sorter}}-{{bam_sort_order}}.bam", **config)
    output:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{{bam_sorter}}-{{bam_sort_order}}.samtools_stats.txt", **config)
    log:
        expand("{log_dir}/samtools_stats/{{sample}}-{{assembly}}-{{bam_sorter}}-{{bam_sort_order}}.log", **config)
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"

def get_featureCounts_bam(wildcards):
    if not 'condition' in samples:
        return expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.samtools-coordinate.bam", **config)
    return expand(f"{{result_dir}}/{{dedup_dir}}/{samples.loc[wildcards.sample, 'condition']}-{wildcards.assembly}.samtools-coordinate.bam", **config)


rule featureCounts:
    """
    Use featureCounts to generate the fraction reads in peaks score (frips/assigned reads).
    https://www.biostars.org/p/337872/
    https://www.biostars.org/p/228636/
    """
    input:
        bam=get_featureCounts_bam,
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
    """
    Generate quality control report for fastq files.
    """
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
                                           "relevant quality control functions."
    qc = []
    for sample in samples[samples['assembly'] == wildcards.assembly].index:
        for function in quality_control:
            qc.extend(function(sample))
    return qc


rule multiqc:
    """
    Aggregate all the quality control metrics for every sample into a single multiqc report.
    """
    input:
       get_qc_files
    output:
        expand("{result_dir}/qc/multiqc_{{assembly}}.html", **config),
        directory(expand("{result_dir}/qc/multiqc_{{assembly}}_data", **config))
    params:
        "{result_dir}/qc/".format(**config)
    log:
        expand("{log_dir}/multiqc_{{assembly}}.log", **config)
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc {input} -o {params} -n multiqc_{wildcards.assembly}.html --config ../../schemas/multiqc_config.yaml > {log} 2>&1"


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
    output = []
    if 'peak_caller' in config:
        if config['peak_caller'] in ['macs2', 'hmmratac']:
            output.append(f"{{result_dir}}/{{dedup_dir}}/{sample}-{samples.loc[sample]['assembly']}.samtools-coordinate.metrics.txt")
        else:
            output.append(f"{{result_dir}}/{{dedup_dir}}/{sample}-{samples.loc[sample]['assembly']}.sambamba-queryname.metrics.txt")
    else:
        output.append(f"{{result_dir}}/{{dedup_dir}}/{sample}-{samples.loc[sample]['assembly']}.{{bam_sorter}}-{{bam_sort_order}}.metrics.txt")

    output.append(f"{{result_dir}}/{{dedup_dir}}/{sample}-{samples.loc[sample]['assembly']}.{{bam_sorter}}-{{bam_sort_order}}.samtools_stats.txt")

    return expand(output, **config)


def get_peak_calling_qc(sample):
    return expand(f"{{result_dir}}/{{peak_caller}}/{sample}-{samples.loc[sample]['assembly']}_featureCounts.txt.summary",
                  **config)
