rule samtools_stats:
    """
    Get general stats from bam files like percentage mapped.
    """
    input:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config)
    output:
        expand("{qc_dir}/dedup/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.samtools_stats.txt", **config)
    log:
        expand("{log_dir}/samtools_stats/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config)
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"


def get_featureCounts_bam(wildcards):
    if wildcards.peak_caller == 'macs2':
        return expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    return expand("{dedup_dir}/{{assembly}}-{{sample}}.sambamba-queryname.bam", **config)


rule featureCounts:
    """
    Use featureCounts to generate the fraction reads in peaks score (frips/assigned reads).
    https://www.biostars.org/p/337872/
    https://www.biostars.org/p/228636/
    """
    input:
        bam=get_featureCounts_bam,
        peak=expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_peaks.narrowPeak", **config)
    output:
        tmp_saf=temp(expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}.saf", **config)),
        real_out=expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_featureCounts.txt", **config),
        summary=expand("{qc_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_featureCounts.txt.summary", **config)
    log:
        expand("{log_dir}/featureCounts/{{assembly}}-{{sample}}-{{peak_caller}}.log", **config)
    threads: 4
    conda:
        "../envs/subread.yaml"
    shell:
        """
        # Make a custom "SAF" file which featureCounts needs:
        awk 'BEGIN{{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' {input.peak} 1> {output.tmp_saf} 2> {log}

        # run featureCounts
        featureCounts -T {threads} -p -a {output.tmp_saf} -F SAF -o {output.real_out} {input.bam} > {log} 2>&1
        
        # move the summary to qc directory
        mv $(dirname {output.real_out})/$(basename {output.summary}) {output.summary}
        """


def get_fastqc_input(wildcards):
    if '_trimmed' in wildcards.fname:
        if 'condition' in samples and config.get('combine_replicates', '') == 'merge' and all(sample not in wildcards.fname for sample in samples.index):
            fqc_input = "{trimmed_dir}/merged/{{fname}}.{fqsuffix}.gz"
        else:
            fqc_input = "{trimmed_dir}/{{fname}}.{fqsuffix}.gz"
    else:
        fqc_input = "{fastq_dir}/{{fname}}.{fqsuffix}.gz"

    return sorted(expand(fqc_input, **config))


rule fastqc:
    """
    Generate quality control report for fastq files.
    """
    input:
        get_fastqc_input
    output:
        f"{config['qc_dir']}/fastqc/{{fname}}_fastqc.html",
        f"{config['qc_dir']}/fastqc/{{fname}}_fastqc.zip"
    log:
        f"{config['log_dir']}/fastqc/{{fname}}.log"
    params:
        f"{config['qc_dir']}/fastqc/"
    conda:
        "../envs/qc.yaml"
    priority: -10
    shell:
        "fastqc {input} -O {params} > {log} 2>&1"


rule InsertSizeMetrics:
    """
    Get the insert size metrics from a (paired-end) bam.
    """
    input:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    output:
        tsv=expand("{qc_dir}/InsertSizeMetrics/{{assembly}}-{{sample}}.tsv", **config),
        pdf=expand("{qc_dir}/InsertSizeMetrics/{{assembly}}-{{sample}}.pdf", **config)
    log:
        f"{config['log_dir']}/InsertSizeMetrics/{{assembly}}-{{sample}}.log"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard CollectInsertSizeMetrics INPUT={input} OUTPUT={output.tsv} H={output.pdf} > {log} 2>&1
        """


def get_qc_files(wildcards):
    assert 'quality_control' in globals(), "When trying to generate multiqc output, make sure that the "\
                                           "variable 'quality_control' exists and contains all the "\
                                           "relevant quality control functions."
    qc = []
    if 'condition' in samples and config.get('combine_replicates', '') == 'merge':
        # trimming qc on individual samples, other qc on merged replicates
        for sample in samples[samples['assembly'] == wildcards.assembly].index:
            qc.extend(get_trimming_qc(sample))

        for condition in set(samples.condition):
            for function in [func for func in quality_control if
                             'get_trimming_qc' is not func.__name__]:
                qc.extend(function(condition))
    else:
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
        expand("{qc_dir}/multiqc_{{assembly}}.html", **config),
        directory(expand("{qc_dir}/multiqc_{{assembly}}_data", **config))
    params:
        dir = "{qc_dir}/".format(**config),
        fqext1 = '_' + config['fqext1']
    log:
        expand("{log_dir}/multiqc_{{assembly}}.log", **config)
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {input} -o {params.dir} -n multiqc_{wildcards.assembly}.html \
        --config ../../schemas/multiqc_config.yaml                           \
        --cl_config "extra_fn_clean_exts: [                                  \
            {{'pattern': ^.*{wildcards.assembly}-, 'type': 'regex'}},        \
            {{'pattern': {params.fqext1},          'type': 'regex'}},        \
            ]" > {log} 2>&1
        """


def get_trimming_qc(sample):
    if config['layout'][sample] == 'SINGLE':
        return expand([f"{{qc_dir}}/fastqc/{sample}_fastqc.zip",
                       f"{{qc_dir}}/fastqc/{sample}_trimmed_fastqc.zip",
                       f"{{qc_dir}}/trimming/{sample}.{{fqsuffix}}.gz_trimming_report.txt"],
                       **config)
    else:
        return expand([f"{{qc_dir}}/fastqc/{sample}_{{fqext}}_fastqc.zip",
                       f"{{qc_dir}}/fastqc/{sample}_{{fqext}}_trimmed_fastqc.zip",
                       f"{{qc_dir}}/trimming/{sample}_{{fqext}}.{{fqsuffix}}.gz_trimming_report.txt"],
                       **config)


def get_alignment_qc(sample):
    output = []
    if 'peak_caller' in config:
        if config['peak_caller'] in ['macs2', 'hmmratac']:
            output.append(f"{{qc_dir}}/dedup/{{{{assembly}}}}-{sample}.samtools-coordinate.metrics.txt")
            output.append(f"{{qc_dir}}/dedup/{{{{assembly}}}}-{sample}.samtools-coordinate.samtools_stats.txt")
        if config['peak_caller'] in ['genrich']:
            output.append(f"{{qc_dir}}/dedup/{{{{assembly}}}}-{sample}.sambamba-queryname.metrics.txt")
            output.append(f"{{qc_dir}}/dedup/{{{{assembly}}}}-{sample}.sambamba-queryname.samtools_stats.txt")
    else:
        output.append(f"{{qc_dir}}/dedup/{{{{assembly}}}}-{sample}.{{bam_sorter}}-{{bam_sort_order}}.metrics.txt")
        output.append(f"{{qc_dir}}/dedup/{{{{assembly}}}}-{sample}.{{bam_sorter}}-{{bam_sort_order}}.samtools_stats.txt")

    if config['layout'][sample] == "PAIRED":
        output.append(f"{{qc_dir}}/InsertSizeMetrics/{{{{assembly}}}}-{sample}.tsv")

    return expand(output, **config)


def get_peak_calling_qc(sample):
    return expand(f"{{qc_dir}}/{{peak_caller}}/{{{{assembly}}}}-{sample}_featureCounts.txt.summary", **config)
