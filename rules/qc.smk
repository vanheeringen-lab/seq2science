import os
import re


rule samtools_stats:
    """
    Get general stats from bam files like percentage mapped.
    """
    input:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-unsieved.bam", **config)
    output:
        expand("{qc_dir}/samtools_stats/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.samtools_stats.txt", **config)
    log:
        expand("{log_dir}/samtools_stats/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config)
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"


# TODO: featurecounts was not set up to accept hmmratac's gappedPeaks,
#  I added this as it was sufficient for rule idr
def get_featureCounts_bam(wildcards):
    if wildcards.peak_caller in ['macs2','hmmratac']:
        return expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    return expand("{dedup_dir}/{{assembly}}-{{sample}}.sambamba-queryname.bam", **config)

# TODO: featurecounts was not set up to accept hmmratac's gappedPeaks,
#  I added this as it was sufficient for rule idr
def get_featureCounts_peak(wildcards):
    peak_sample = brep_from_trep[wildcards.sample]
    ftype = get_ftype(wildcards.peak_caller)
    return expand(f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{peak_sample}_peaks.{ftype}", **config)

rule featureCounts:
    """
    Use featureCounts to generate the fraction reads in peaks score (frips/assigned reads).
    https://www.biostars.org/p/337872/
    https://www.biostars.org/p/228636/
    """
    input:
        bam=get_featureCounts_bam,
        peak=get_featureCounts_peak,
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
        if 'replicate' in samples and config.get('technical_replicates', '') == 'merge' and all(sample not in wildcards.fname for sample in samples.index):
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

def get_chrM_name(wildcards, input):
    if os.path.exists(str(input.chr_names)):
        name = [chrm for chrm in ['chrM', 'MT'] if chrm in open(str(input.chr_names), 'r').read()]
        if len(name) > 0:
            return name[0]
    return "no_chrm_found"


rule MTNucRatioCalculator:
    """
    Calculate the ratio mitochondrial dna in your sample.
    """
    input:
        bam=expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-unsieved.bam", **config),
        chr_names=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json", **config)
    log:
        f"{config['log_dir']}/MTNucRatioCalculator/{{assembly}}-{{sample}}.log"
    conda:
        "../envs/mtnucratio.yaml"
    params:
        lambda wildcards, input: get_chrM_name(wildcards, input)
    shell:
        """
        mtnucratio {input.bam} {params}
        """


def fingerprint_multiBamSummary_input(wildcards):
    output = {"bams": set(), "bais": set()}

    for trep in set(treps[treps['assembly'] == wildcards.assembly].index):
        output["bams"].update(expand(f"{{dedup_dir}}/{wildcards.assembly}-{trep}.samtools-coordinate.bam", **config))
        output["bais"].update(expand(f"{{dedup_dir}}/{wildcards.assembly}-{trep}.samtools-coordinate.bam.bai", **config))
        if "control" in treps and isinstance(treps.loc[trep, "control"], str):
            control = treps.loc[trep, "control"]
            output["bams"].update(expand(f"{{dedup_dir}}/{wildcards.assembly}-{control}.samtools-coordinate.bam", **config))
            output["bais"].update(expand(f"{{dedup_dir}}/{wildcards.assembly}-{control}.samtools-coordinate.bam.bai", **config))

    return output


def get_descriptive_names(wildcards, input):
    if "descriptive_name" not in treps:
        return ""

    labels = ""
    for file in input:
        trep = re.findall(f"{wildcards.assembly}-([^\.]+)\.+?", file.split("/")[-1])[0]
        if "control" in treps and trep not in treps.index:
            labels += f"control_{trep} "
        else:
            labels += treps.loc[trep, "descriptive_name"] + " "

    return labels


rule plotFingerprint:
    input:
        unpack(fingerprint_multiBamSummary_input)
    output:
        expand("{qc_dir}/plotFingerprint/{{assembly}}.tsv", **config)
    log:
        expand("{log_dir}/plotFingerprint/{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/plotFingerprint/{{{{assembly}}}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 16
    params:
        lambda wildcards, input: "--labels " + get_descriptive_names(wildcards, input.bams) if
                                 get_descriptive_names(wildcards, input.bams) != "" else ""
    shell:
        """
        plotFingerprint -b {input.bams} {params} --outRawCounts {output} -p {threads} > {log} 2>&1 
        """


def computematrix_input(wildcards):
    output = []

    for trep in set(treps[treps['assembly'] == assembly].index):
        output.append(expand(f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{trep}.bw", **config)[0])

    return output


rule computeMatrix:
    input:
        bw=computematrix_input
    output:
        expand("{qc_dir}/computeMatrix/{{assembly}}-{{peak_caller}}.mat.gz", **config)
    log:
        expand("{log_dir}/computeMatrix/{{assembly}}-{{peak_caller}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/computeMatrix/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 16
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    params:
        labels=lambda wildcards, input: "--samplesLabel " + get_descriptive_names(wildcards, input.bw) if
                                 get_descriptive_names(wildcards, input.bw) != "" else "",
        annotation=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config)  # TODO: move genomepy to checkpoint and this as input
    shell:
        """
        computeMatrix scale-regions -S {input.bw} {params.labels} -R {params.annotation} -p {threads} -b 2000 -a 500 -o {output} > {log} 2>&1
        """


rule plotProfile:
    input:
        rules.computeMatrix.output
    output:
        file=expand("{qc_dir}/plotProfile/{{assembly}}-{{peak_caller}}.tsv", **config),
        img=expand("{qc_dir}/plotProfile/{{assembly}}-{{peak_caller}}.png", **config)
    log:
        expand("{log_dir}/plotProfile/{{assembly}}-{{peak_caller}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/plotProfile/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    shell:
        """
        plotProfile -m {input} --outFileName {output.img} --outFileNameData {output.file} > {log} 2>&1
        """


rule multiBamSummary:
    input:
        unpack(fingerprint_multiBamSummary_input)
    output:
        expand("{qc_dir}/multiBamSummary/{{assembly}}.npz", **config)
    log:
        expand("{log_dir}/multiBamSummary/{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/multiBamSummary/{{assembly}}.benchmark.txt", **config)[0]
    threads: 16
    params:
        lambda wildcards, input: "--labels " + get_descriptive_names(wildcards, input.bams) if
                                 get_descriptive_names(wildcards, input.bams) != "" else "",
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    shell:
        """
        multiBamSummary bins --bamfiles {input.bams} -out {output} {params} -p {threads} > {log} 2>&1
        """


rule plotCorrelation:
    input:
        rules.multiBamSummary.output
    output:
        expand("{qc_dir}/plotCorrelation/{{assembly}}.tsv", **config)
    log:
        expand("{log_dir}/plotCorrelation/{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/plotCorrelation/{{assembly}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    shell:
        """
        plotCorrelation --corData {input} --outFileCorMatrix {output} -c spearman -p heatmap > {log} 2>&1
        """


rule plotPCA:
    input:
        rules.multiBamSummary.output
    output:
        expand("{qc_dir}/plotPCA/{{assembly}}.tsv", **config)
    log:
        expand("{log_dir}/plotPCA/{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/plotPCA/{{assembly}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    shell:
        """
        plotPCA --corData {input} --outFileNameData {output} > {log} 2>&1
        """


rule multiqc_header_info:
    """
    Generate a multiqc header file with contact info and date of multiqc generation.
    """
    output:
        temp(expand('{qc_dir}/header_info.yaml', **config))
    run:
        import os
        import copy
        from datetime import date

        cwd = os.getcwd().split('/')[-1]
        mail = config.get('email', 'none@provided.com')
        date = date.today().strftime("%B %d, %Y")

        with open(output[0], "w") as f:
            f.write(f"report_header_info:\n"
                    f"    - Contact E-mail: '{mail}'\n"
                    f"    - Workflow: '{cwd}'\n"
                    f"    - Date: '{date}'\n")


rule multiqc_rename_buttons:
    """
    Generate rename buttons.
    """
    output:
        temp(expand('{qc_dir}/sample_names_{{assembly}}.tsv', **config))
    run:
        newsamples = samples[samples["assembly"] == wildcards.assembly].reset_index(level=0, inplace=False)
        newsamples = newsamples.drop(["assembly"], axis=1)
        newsamples.to_csv(output[0], sep="\t", index=False)


def get_qc_files(wildcards):
    assert 'quality_control' in globals(), "When trying to generate multiqc output, make sure that the "\
                                           "variable 'quality_control' exists and contains all the "\
                                           "relevant quality control functions."
    qc = dict()
    qc['header'] = expand('{qc_dir}/header_info.yaml', **config)[0]
    qc['sample_names'] = expand('{qc_dir}/sample_names_{{assembly}}.tsv', **config)[0]
    qc['files'] = set()

    # trimming qc on individual samples
    if get_trimming_qc in quality_control:
        for sample in samples[samples['assembly'] == wildcards.assembly].index:
            qc['files'].update(get_trimming_qc(sample))

    # qc on merged technical replicates/samples
    if get_alignment_qc in quality_control:
        for replicate in treps[treps['assembly'] == wildcards.assembly].index:
            for function in [func for func in quality_control if
                             func.__name__ not in ['get_peak_calling_qc', 'get_trimming_qc']]:
                qc['files'].update(function(replicate))

    # qc on combined biological replicates/samples
    if get_peak_calling_qc in quality_control:
        for trep in treps[treps['assembly'] == wildcards.assembly].index:
            qc['files'].update(get_peak_calling_qc(trep))

    return qc


rule multiqc:
    """
    Aggregate all the quality control metrics for every sample into a single multiqc report.
    """
    input:
        unpack(get_qc_files)
    output:
        expand("{qc_dir}/multiqc_{{assembly}}.html", **config),
        directory(expand("{qc_dir}/multiqc_{{assembly}}_data", **config))
    params:
        dir = "{qc_dir}/".format(**config),
        fqext1 = '_' + config['fqext1'],
    log:
        expand("{log_dir}/multiqc_{{assembly}}.log", **config)
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {input.files} -o {params.dir} -n multiqc_{wildcards.assembly}.html \
        --config ../../schemas/multiqc_config.yaml                                 \
        --config {input.header}                                                    \
        --sample-names {input.sample_names}                                        \
        --cl_config "extra_fn_clean_exts: [                                        \
            {{'pattern': ^.*{wildcards.assembly}-, 'type': 'regex'}},              \
            {{'pattern': {params.fqext1},          'type': 'regex'}},              \
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
    output.append(f"{{qc_dir}}/samtools_stats/{{{{assembly}}}}-{sample}.samtools-coordinate.samtools_stats.txt")

    if "atac_seq" in get_workflow():
        # add insert size metrics
        if config['layout'][sample] == "PAIRED":
            output.append(f"{{qc_dir}}/InsertSizeMetrics/{{{{assembly}}}}-{sample}.tsv")

    # get the ratio mitochondrial dna
    output.append(f"{{result_dir}}/{config['aligner']}/{{{{assembly}}}}-{sample}.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json")

    if get_workflow() in ["alignment", "chip_seq", "atac_seq"]:
        output.append("{qc_dir}/plotFingerprint/{{assembly}}.tsv")
    if len(breps["assembly"] == treps.loc[sample, "assembly"]) > 1:
        output.append("{qc_dir}/plotCorrelation/{{assembly}}.tsv")
        output.append("{qc_dir}/plotPCA/{{assembly}}.tsv")

    return expand(output, **config)


def get_peak_calling_qc(sample):
    output = []

    # add frips score (featurecounts)
    output.extend(expand(f"{{qc_dir}}/{{peak_caller}}/{{{{assembly}}}}-{sample}_featureCounts.txt.summary", **config))

    # macs specific QC
    if "macs2" in config["peak_caller"]:
        output.extend(expand(f"{{result_dir}}/macs2/{{{{assembly}}}}-{sample}_peaks.xls", **config))

    # deeptools profile
    assembly = treps.loc[sample, "assembly"]
    # TODO: replace with genomepy checkpoint in the future
    if assembly.lower() in ["ce10", "dm3", "hg38", "hg19", "mm9", "mm10", "grch38.p13", "grch38"]:
        output.extend(expand("{qc_dir}/plotProfile/{{assembly}}-{peak_caller}.tsv", **config))

    return output

