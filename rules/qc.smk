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


rule MTNucRatioCalculator:
    """
    Calculate the ratio mitochondrial dna in your sample.
    Note: version 0.5 is bugged in what output it uses, and currently the {output} is just there
    for form.
    """
    input:
        bam=expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-unsieved.bam", **config),
        chr_names=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.bam.mtnucratiomtnuc.json", **config)
    conda:
        "../envs/mtnucratio.yaml"
    params:
        lambda wildcards, input: [chrm for chrm in ['chrM', 'MT'] if chrm in open(str(input.chr_names), 'r').read()][0]
    shell:
        """
        mtnucratio {input.bam} {output} {params}
        """


rule multiqc_header_info:
    """
    Generate a multiqc header file with contact info and date of multiqc generation.
    """
    output:
        expand('{qc_dir}/header_info.yaml', **config)
    run:
        import os
        import copy
        from datetime import date

        cwd = os.getcwd().split('/')[-1]
        mail = config.get('email', 'none@provided.com')
        date = date.today().strftime("%B %d, %Y")

        description = ['        <dl class=dl-horizontal>']
        for key, value in config.items():
            if not any([i in key for i in ['dir', 'ascp', 'ncbi']]):
                # description.append(f'{key}:\t\t{value}<br /> ')
                description.append(f'            <dt>{key}</dt><dd>{value}</dd>')
        description.append('        </dl>')

        # make a table of our samples
        samples_table = copy.copy(samples)
        samples_table.index.name = None
        samples_table = samples_table.to_html().split('\n')
        samples_table = ['        ' + sample for sample in samples_table]
        description.extend(samples_table)

        description = '\n'.join(description)
        description = description.replace('\'', '').replace(':', '\:')

        with open(output[0], "w") as f:
            f.write(f"report_header_info:\n"
                    f"    - Contact E-mail: '{mail}'\n"
                    f"    - Workflow: '{cwd}'\n"
                    f"    - Date: '{date}'\n"
                    f"\n\n"
                    f"custom_data:\n"
                    f"    id: 'Whatup fellas'\n"
                    f"    section_name: 'Let me tell you a story'\n"
                    f"    plot_type: 'html'\n"
                    f"    data: |\n"
                    f"        The pipeline was run with these configuration options and samples\:\n"
                    f"{description}\n"
                    # f"custom_data:\n"
                    # f"    my_data_type:\n"
                    # f"        id: 'mqc_config_file_section'\n"
                    # f"        section_name: 'Configuration options<br /> '\n"
                    # f"        description: {description}\n"
                    # f"        plot_type: 'bargraph'\n"
                    # f"        data:\n"
                    # f"            sample_a:\n"
                    # f"                empty: 100\n"
                    )


rule multiqc_config_info:
    """
    Generate a multiqc header file with contact info and date of multiqc generation.
    """
    output:
        expand('{qc_dir}/config_info.yaml', **config)
    run:
        with open(output[0], "w") as f:
            f.write(
                "# section_name: 'Custom data file'"
                "# description: 'This output is described in the file header. Any MultiQC installation will understand it without prior configuration.'"
                "# format: 'csv'"
                "# plot_type: 'table'"
                "# pconfig:"
                "#    id: 'custom_bargraph_w_header'"
                "#    ylab: 'Number of things'"
                "Sample Name,some feature"
                "GSM2219688_pass_1,42")


def get_qc_files(wildcards):
    assert 'quality_control' in globals(), "When trying to generate multiqc output, make sure that the "\
                                           "variable 'quality_control' exists and contains all the "\
                                           "relevant quality control functions."
    qc = []
    if 'condition' in samples and config.get('combine_replicates', '') == 'merge':
        # trimming qc on individual samples, other qc on merged replicates
        for sample in samples[samples['assembly'] == wildcards.assembly].index:
            qc.extend(get_trimming_qc(sample))

        for condition in set(samples[samples['assembly'] == wildcards.assembly].condition):
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
        qc=get_qc_files,
        header=expand('{qc_dir}/header_info.yaml', **config),
        # config=expand('{qc_dir}/config_info.yaml', **config)
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
        --config {input.header}                                              \
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

    # add samtools stats
    output.append(f"{{qc_dir}}/dedup/{{{{assembly}}}}-{sample}.samtools-coordinate.metrics.txt")
    output.append(f"{{qc_dir}}/dedup/{{{{assembly}}}}-{sample}.samtools-coordinate.samtools_stats.txt")

    # add insert size metrics
    if config['layout'][sample] == "PAIRED":
        output.append(f"{{qc_dir}}/InsertSizeMetrics/{{{{assembly}}}}-{sample}.tsv")

    # get the ratio mitochondrial dna
    assembly = samples.loc[sample, 'assembly']
    checkpoints.get_genome.get(assembly=assembly)
    if len([chrm for chrm in ['chrM', 'MT'] if chrm in open(f"{config['genome_dir']}/{assembly}/{assembly}.fa.sizes", 'r').read()]) > 0:
        output.append(f"{{result_dir}}/{config['aligner']}/{{{{assembly}}}}-{sample}.samtools-coordinate.bam.mtnucratiomtnuc.json")
    return expand(output, **config)


def get_peak_calling_qc(sample):
    return expand(f"{{qc_dir}}/{{peak_caller}}/{{{{assembly}}}}-{sample}_featureCounts.txt.summary", **config)
