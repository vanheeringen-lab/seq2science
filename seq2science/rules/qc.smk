import re

import seq2science
from seq2science.util import sieve_bam, get_bustools_rid


localrules: multiqc_header_info, multiqc_rename_buttons, multiqc_filter_buttons, multiqc_samplesconfig, multiqc_schema, combine_qc_files

def samtools_stats_input(wildcards):
    if wildcards.directory == config["aligner"]:
        return expand("{result_dir}/{{directory}}/{{assembly}}-{{sample}}.samtools-coordinate-unsieved.bam", **config)
    return expand("{final_bam_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config)


rule samtools_stats:
    """
    Get general stats from bam files like percentage mapped.
    """
    input:
        samtools_stats_input
    output:
        expand("{qc_dir}/samtools_stats/{{directory}}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.samtools_stats.txt", **config)
    log:
        expand("{log_dir}/samtools_stats/{{directory}}/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config)
    message: explain_rule("samtools_stats")
    resources:
        time="0-06:00:00"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"


# TODO: featurecounts was not set up to accept hmmratac's gappedPeaks,
#  I added this as it was sufficient for rule idr
def get_featureCounts_bam(wildcards):
    if wildcards.peak_caller in ['macs2','hmmratac']:
        return expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    return expand("{final_bam_dir}/{{assembly}}-{{sample}}.sambamba-queryname.bam", **config)

# TODO: featurecounts was not set up to accept hmmratac's gappedPeaks,
#  I added this as it was sufficient for rule idr
def get_featureCounts_peak(wildcards):
    peak_sample = brep_from_trep[wildcards.sample]
    ftype = get_peak_ftype(wildcards.peak_caller)
    return expand(f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{peak_sample}_peaks.{ftype}", **config)

rule featureCounts:
    """
    Use featureCounts to generate the fraction reads in peaks score 
    (frips/assigned reads). See: https://www.biostars.org/p/337872/ and
    https://www.biostars.org/p/228636/
    """
    input:
        bam=get_featureCounts_bam,
        peak=get_featureCounts_peak,
    output:
        tmp_saf=temp(expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}.saf", **config)),
        real_out=expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_featureCounts.txt", **config),
        summary=expand("{qc_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_featureCounts.txt.summary", **config)
    message: explain_rule("featureCounts_qc")
    log:
        expand("{log_dir}/featureCounts/{{assembly}}-{{sample}}-{{peak_caller}}.log", **config)
    threads: 4
    conda:
        "../envs/subread.yaml"
    shell:
        """
        # Make a custom "SAF" file which featureCounts needs:
        awk 'BEGIN{{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' \
        {input.peak} 1> {output.tmp_saf} 2> {log}

        # run featureCounts
        featureCounts -T {threads} -p -a {output.tmp_saf} -F SAF -o {output.real_out} {input.bam} > {log} 2>&1
        
        # move the summary to qc directory
        mv $(dirname {output.real_out})/$(basename {output.summary}) {output.summary}
        """


def get_fastqc_input(wildcards):
    if '_trimmed' in wildcards.fname:
        fqc_input = "{trimmed_dir}/{{fname}}.{fqsuffix}.gz"
    else:
        fqc_input = "{fastq_dir}/{{fname}}.{fqsuffix}.gz"

    return sorted(expand(fqc_input, **config))


if config["trimmer"] == "trimgalore":
    rule fastqc:
        """
        Generate quality control report for fastq files.
        """
        input:
            get_fastqc_input
        output:
            f"{config['qc_dir']}/fastqc/{{fname}}_fastqc.html",
            f"{config['qc_dir']}/fastqc/{{fname}}_fastqc.zip"
        message:explain_rule("fastqc")
        log:
            f"{config['log_dir']}/fastqc/{{fname}}.log"
        params:
            f"{config['qc_dir']}/fastqc/"
        conda:
            "../envs/fastqc.yaml"
        priority: -10
        shell:
            """
            fastqc {input} -O {params} > {log} 2>&1
            """


elif config["trimmer"] == "fastp" or "kmer_estimation" in config:
    ruleorder: fastp_qc_PE> fastp_qc_SE > fastp_PE > fastp_SE


    checkpoint fastp_qc_SE:
        """
        Get quality scores for (technical) replicates
        """
        input:
            expand("{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config),
        output:
            qc_json=expand("{qc_dir}/trimming/{{sample}}.fastp.json", **config),
            qc_html=expand("{qc_dir}/trimming/{{sample}}.fastp.html", **config),
        conda:
            "../envs/fastp.yaml"
        threads: 1
        wildcard_constraints:
            sample="|".join(merged_treps_single) if len(merged_treps_single) else "$a"
        log:
            expand("{log_dir}/fastp_qc_SE/{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/fastp_qc_SE/{{sample}}.benchmark.txt", **config)[0]
        priority: -10
        params:
            fqsuffix=config["fqsuffix"],
        shell:
            """\
            fastp -w {threads} --in1 {input} -h {output.qc_html} -j {output.qc_json} > {log} 2>&1
            """


    rule fastp_qc_PE:
        """
        Get quality scores for (technical) replicates
        """
        input:
            r1=expand("{trimmed_dir}/{{sample}}_{fqext1}_trimmed.{fqsuffix}.gz", **config),
            r2=expand("{trimmed_dir}/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config),
        output:
            qc_json=expand("{qc_dir}/trimming/{{sample}}.fastp.json", **config),
            qc_html=expand("{qc_dir}/trimming/{{sample}}.fastp.html", **config),
        conda:
            "../envs/fastp.yaml"
        threads: 1
        wildcard_constraints:
            sample="|".join(merged_treps_paired) if len(merged_treps_paired) else "$a"
        log:
            expand("{log_dir}/fastp_qc_PE/{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/fastp_qc_PE/{{sample}}.benchmark.txt", **config)[0]
        priority: -10
        params:
            config=config["trimoptions"],
        shell:
            """\
            fastp -w {threads} --in1 {input[0]} --in2 {input[1]} \
            -h {output.qc_html} -j {output.qc_json} \
            {params.config} > {log} 2>&1
            """


rule insert_size_metrics:
    """
    Get insert size metrics from a (paired-end) bam. This score is then used by
    MultiQC in the report.
    """
    input:
        expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    output:
        tsv=expand("{qc_dir}/InsertSizeMetrics/{{assembly}}-{{sample}}.tsv", **config),
        pdf=expand("{qc_dir}/InsertSizeMetrics/{{assembly}}-{{sample}}.pdf", **config)
    log:
        f"{config['log_dir']}/InsertSizeMetrics/{{assembly}}-{{sample}}.log"
    conda:
        "../envs/picard.yaml"
    wildcard_constraints:
        sample=".+",
    resources:
        time="0-06:00:00"
    shell:
        """
        picard CollectInsertSizeMetrics INPUT={input} \
        OUTPUT={output.tsv} H={output.pdf} > {log} 2>&1
        """

def get_chrM_name(wildcards, input):
    if os.path.exists(str(input.chr_names)):
        name = [chrm for chrm in ['chrM', 'MT'] if chrm in open(str(input.chr_names), 'r').read()]
        if len(name) > 0:
            return name[0]
    return "no_chrm_found"


rule mt_nuc_ratio_calculator:
    """
    Estimate the amount of nuclear and mitochondrial reads in a sample. This metric
    is especially important in ATAC-seq experiments where mitochondrial DNA can
    be overrepresented. Reads are classified as mitochondrial reads if they map
    against either "chrM" or "MT".
    
    These values are aggregated and displayed in the MultiQC report.
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
        mitochondria=lambda wildcards, input: get_chrM_name(wildcards, input)
    resources:
        time="0-06:00:00"
    shell:
        """
        mtnucratio {input.bam} {params.mitochondria}
        """


def fingerprint_multiBamSummary_input(wildcards):
    output = {"bams": set(), "bais": set()}

    for trep in set(treps[treps['assembly'] == ori_assembly(wildcards.assembly)].index):
        output["bams"].update(expand(f"{{final_bam_dir}}/{wildcards.assembly}-{trep}.samtools-coordinate.bam", **config))
        output["bais"].update(expand(f"{{final_bam_dir}}/{wildcards.assembly}-{trep}.samtools-coordinate.bam.bai", **config))
        if "control" in treps and isinstance(treps.loc[trep, "control"], str):
            control = treps.loc[trep, "control"]
            output["bams"].update(expand(f"{{final_bam_dir}}/{wildcards.assembly}-{control}.samtools-coordinate.bam", **config))
            output["bais"].update(expand(f"{{final_bam_dir}}/{wildcards.assembly}-{control}.samtools-coordinate.bam.bai", **config))

    return output


def get_descriptive_names(wildcards, input):
    if "descriptive_name" not in samples:
        return ""

    labels = ""
    for file in input:
        # extract trep name from filepath (between assembly- and .sorter)
        trep = file[file.rfind(f"{wildcards.assembly}-"):].replace(f"{wildcards.assembly}-", "")

        # get the sample name
        if trep.find(".sam") != -1:
            trep = trep[:trep.find(".sam")]
        elif trep.find(".bw") != -1:
            trep = trep[:trep.find(".bw")]
        elif trep.find("_summits.bed") != -1:
            trep = trep[:trep.find("_summits.bed")]
        else:
            raise ValueError

        if trep in breps.index:
            if len(treps_from_brep[(trep, wildcards.assembly)]) == 1 and trep in samples.index:
                labels += samples.loc[trep, "descriptive_name"] + " "
            else:
                labels += trep + " "
        elif "control" in treps and trep not in treps.index:
            labels += f"control_{trep} "
        elif trep in samples.index:
            labels += samples.loc[trep, "descriptive_name"] + " "
        else:
            labels += trep + " "

    return labels


rule plotFingerprint:
    """
    Plot the "fingerprint" of your bams, using deeptools. 
    """
    input:
        unpack(fingerprint_multiBamSummary_input)
    output:
        expand("{qc_dir}/plotFingerprint/{{assembly}}.tsv", **config)
    log:
        expand("{log_dir}/plotFingerprint/{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/plotFingerprint/{{assembly}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 16
    resources:
        mem_gb=5,
    params:
        lambda wildcards, input: "--labels " + get_descriptive_names(wildcards, input.bams) if
                                 get_descriptive_names(wildcards, input.bams) != "" else ""
    shell:
        """
        plotFingerprint -b {input.bams} {params} --outRawCounts {output} -p {threads} > {log} 2>&1 
        """


def computeMatrix_input(wildcards):
    output = []

    treps_seen = set()
    for trep in treps[treps['assembly'] == ori_assembly(wildcards.assembly)].index:
        if trep not in treps_seen:
            output.append(expand(f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{trep}.bw", **config)[0])
        treps_seen.add(trep)

    return output


rule computeMatrix_gene:
    """
    Pre-compute correlations between bams using deeptools.
    """
    input:
        bw=computeMatrix_input
    output:
        expand("{qc_dir}/computeMatrix_gene/{{assembly}}-{{peak_caller}}.mat.gz", **config)
    log:
        expand("{log_dir}/computeMatrix_gene/{{assembly}}-{{peak_caller}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/computeMatrix_gene/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 16
    resources:
        deeptools_limit=lambda wildcards, threads: threads,
        mem_gb=50,
    params:
        labels=lambda wildcards, input: "--samplesLabel " + get_descriptive_names(wildcards, input.bw) if
                                 get_descriptive_names(wildcards, input.bw) != "" else "",
        annotation=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config)  # TODO: move genomepy to checkpoint and this as input
    shell:
        """
        computeMatrix scale-regions -S {input.bw} {params.labels} -R {params.annotation} \
        -p {threads} -b 2000 -a 500 -o {output} > {log} 2>&1
        """


rule plotProfile_gene:
    """
    Plot the gene profile using deeptools.
    """
    input:
        rules.computeMatrix_gene.output
    output:
        file=expand("{qc_dir}/plotProfile_gene/{{assembly}}-{{peak_caller}}.tsv", **config),
        img=expand("{qc_dir}/plotProfile_gene/{{assembly}}-{{peak_caller}}.png", **config)
    log:
        expand("{log_dir}/plotProfile_gene/{{assembly}}-{{peak_caller}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/plotProfile_gene/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    shell:
        """
        plotProfile -m {input} --outFileName {output.img} --outFileNameData {output.file} > {log} 2>&1
        """


rule computeMatrix_peak:
    """
    Pre-compute correlations between bams using deeptools.
    """
    input:
        bw=computeMatrix_input,
        peaks=expand("{qc_dir}/computeMatrix_peak/{{assembly}}-{{peak_caller}}_N{{nrpeaks}}.bed", **config)
    output:
        expand("{qc_dir}/computeMatrix_peak/{{assembly}}-{{peak_caller}}_N{{nrpeaks}}.mat.gz", **config)
    log:
        expand("{log_dir}/computeMatrix_peak/{{assembly}}-{{peak_caller}}_N{{nrpeaks}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/computeMatrix_peak/{{assembly}}-{{peak_caller}}_N{{nrpeaks}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 16
    resources:
        deeptools_limit=lambda wildcards, threads: threads,
        mem_gb=50,
    params:
        labels=lambda wildcards, input: "--samplesLabel " + get_descriptive_names(wildcards, input.bw) if
                                 get_descriptive_names(wildcards, input.bw) != "" else "",
    shell:
        """
        computeMatrix scale-regions -S {input.bw} {params.labels} -R {input.peaks} \
        -p {threads} --binSize 5 -o {output} > {log} 2>&1
        """


rule plotHeatmap_peak:
    """
    Plot the peak heatmap using deepTools.
    """
    input:
        rules.computeMatrix_peak.output
    output:
        img=expand("{qc_dir}/plotHeatmap_peaks/N{{nrpeaks}}-{{assembly}}-deepTools_{{peak_caller}}_heatmap_mqc.png", **config),
    log:
        expand("{log_dir}/plotHeatmap_peaks/{{assembly}}-{{peak_caller}}_N{{nrpeaks}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/plotHeatmap_peaks/{{assembly}}-{{peak_caller}}_N{{nrpeaks}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    params: 
        params=config.get("heatmap_deeptools_options", ""),
        slop=config.get("heatmap_slop", 0)
    shell:
        """
        plotHeatmap --matrixFile {input} --outFileName {output.img} {params.params} --startLabel="-{params.slop}" --endLabel={params.slop} > {log} 2>&1
        """


rule multiBamSummary:
    """
    Pre-compute a bam summary with deeptools.
    """
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
        names=lambda wildcards, input: "--labels " + get_descriptive_names(wildcards, input.bams) if
                                       get_descriptive_names(wildcards, input.bams) != "" else "",
        params=config["deeptools_multibamsummary"]
    message: explain_rule("computeMatrix")
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads,
        mem_gb=4,
    shell:
        """
        multiBamSummary bins --bamfiles {input.bams} -out {output} {params.names} \
        {params.params} -p {threads} > {log} 2>&1
        """


rule plotCorrelation:
    """
    Calculate the correlation between bams with deepTools.
    """
    input:
        rules.multiBamSummary.output
    output:
        expand("{qc_dir}/plotCorrelation/{{assembly}}-deepTools_{{method}}_correlation_clustering_mqc.png", **config)
    log:
        expand("{log_dir}/plotCorrelation/{{assembly}}-{{method}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/plotCorrelation/{{assembly}}-{{method}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads
    params:
        title=lambda wildcards: ('"Spearman Correlation of Read Counts"'
                                 if wildcards.method == "spearman" else
                                 '"Pearson Correlation of Read Counts"'),
        params=config["deeptools_plotcorrelation"]
    shell:
        """
        plotCorrelation --corData {input} --plotFile {output} -c {wildcards.method} \
        -p heatmap --plotTitle {params.title} {params.params} > {log} 2>&1
        """


rule plotPCA:
    """
    Plot a PCA between bams using deepTools.
    """
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

def get_summits_bed(wildcards):
    return expand(
        [
            f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_summits.bed"
            for replicate in breps[breps["assembly"] == ori_assembly(wildcards.assembly)].index
        ],
        **config,
    )


rule chipseeker:
    input:
        narrowpeaks=get_summits_bed,
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config)
    output:
        img1=expand("{qc_dir}/chipseeker/{{assembly}}-{{peak_caller}}_img1_mqc.png", **config),
        img2=expand("{qc_dir}/chipseeker/{{assembly}}-{{peak_caller}}_img2_mqc.png", **config),
    params:
        names=lambda wildcards, input: get_descriptive_names(wildcards, input.narrowpeaks)
    log:
        expand("{log_dir}/chipseeker/{{assembly}}-{{peak_caller}}.log", **config)
    conda:
        "../envs/chipseeker.yaml"
    message: explain_rule("chipseeker")
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/chipseeker.R"


rule multiqc_header_info:
    """
    Generate a multiqc header file with contact info and date of multiqc generation.
    """
    output:
        temp(expand('{qc_dir}/header_info.yaml', **config))
    params:
        config={k: v for k, v in config.items() if k not in ["no_config_log", "cli_call"]}  # helps resolve changed config options, ignore no_config_log
    run:
        import os
        from datetime import date

        workflow = get_workflow().replace("_", "-")
        date = date.today().strftime("%B %d, %Y")
        project = os.getcwd().split('/')[-1]
        mail = config.get('email', 'none@provided.com')

        with open(output[0], "w") as f:
            f.write(f"report_header_info:\n"
                    f"    - Workflow: '{workflow}'\n"
                    f"    - Date: '{date}'\n"
                    f"    - Project: '{project}'\n"
                    f"    - Contact E-mail: '{mail}'\n"
            )


rule multiqc_explain:
    output:
        expand('{log_dir}/workflow_explanation_mqc.html', **config)
    params:
        config={k: v for k, v in config.items() if k not in ["no_config_log", "cli_call"]}  # helps resolve changed config options, ignore no_config_log
    run:
        import subprocess

        # make a copy of the cli_call
        cli_call = config["cli_call"][:]

        # convert run into explain
        cli_call[1] = "explain"

        # remove run specific commands
        cli_call = [arg for arg in cli_call if arg not in {"--unlock", "--rerun-incomplete", "-k", "--keep-going", "--skip-rerun", "--reason", "-r", "--dryrun", "-n"}]
        cores_index = cli_call.index("--cores" if "--cores" in cli_call else "-j")
        del cli_call[cores_index:cores_index+2]
        
        # make sure to get pretty hyperrefs
        cli_call.append("--hyperref")

        result = subprocess.run(cli_call, stdout=subprocess.PIPE)
        explanation = result.stdout.decode('utf-8')

        with open(output[0], "w") as f:
            f.write("<!--\n" 
                    "id: 'explanation'\n" 
                    "section_name: 'Workflow explanation'\n" 
                    "-->\n"
                    f"{explanation}")


rule multiqc_rename_buttons:
    """
    Generate rename buttons for the multiqc report.
    """
    output:
        temp(expand('{qc_dir}/sample_names_{{assembly}}.tsv', **config))
    params:
        samples  # helps resolve changed params if e.g. descriptive names change
    run:
        newsamples = samples[samples["assembly"] == ori_assembly(wildcards.assembly)].reset_index(level=0, inplace=False)
        newsamples = newsamples.drop(["assembly"], axis=1)
        newsamples.to_csv(output[0], sep="\t", index=False)


rule multiqc_filter_buttons:
    """
    Generate filter buttons.
    """
    output:
        temp(expand('{qc_dir}/sample_filters_{{assembly}}.tsv', **config))
    params:
        samples  # helps resolve changed params if e.g. descriptive names change
    run:
        with open(output[0], "w") as f:
            f.write("Read Group 1 & Alignment\thide\t_R2\n"
                    "Read Group 2\tshow\t_R2\n")


rule multiqc_samplesconfig:
    """
    Add a section in the multiqc report that reports the samples.tsv and config.yaml
    """
    output:
        temp(expand('{qc_dir}/samplesconfig_mqc.html', **config))
    params:
        config_used=len(workflow.overwrite_configfiles) > 0,
        configfile=workflow.overwrite_configfiles[-1],
        sanitized_samples=sanitized_samples,  # helps resolve changed config options
        config={k: v for k, v in config.items() if k not in ["no_config_log", "cli_call"]}  # helps resolve changed config options, ignore no_config_log
    conda:
        "../envs/htmltable.yaml"
    script:
        f"{config['rule_dir']}/../scripts/multiqc_samplesconfig.py"

rule multiqc_schema:
    """
    Generate a multiqc config schema. Used for the ordering of modules.
    """
    output:
        temp(expand('{qc_dir}/schema.yaml', **config))
    run:
        with open(f'{config["rule_dir"]}/../schemas/multiqc_config.yaml') as cookie_schema:
            cookie = cookie_schema.read()
        
        while len(list(re.finditer("{{.*}}", cookie))):
            match_obj = next(re.finditer("{{.*}}", cookie))
            span = match_obj.span()
            match = eval(match_obj.group(0)[2:-2])  # without the curly braces
            cookie = cookie[:span[0]] + match + cookie[span[1]:]

        with open(output[0], "w") as out_file:
            out_file.write(cookie)


def get_qc_files(wildcards):
    qc = dict()
    assert 'quality_control' in globals(), "When trying to generate multiqc output, make sure that the "\
                                           "variable 'quality_control' exists and contains all the "\
                                           "relevant quality control functions."
    qc['files'] = set(expand(['{qc_dir}/samplesconfig_mqc.html',
                              '{log_dir}/workflow_explanation_mqc.html'], **config))

    # trimming qc on individual samples
    if get_trimming_qc in quality_control:
        # scatac seq only on treps, not on single samples
        if get_workflow() != "scatac_seq":
            for sample in samples[samples['assembly'] == ori_assembly(wildcards.assembly)].index:
                qc['files'].update(get_trimming_qc(sample))

    # qc on merged technical replicates/samples
    if get_alignment_qc in quality_control:
        for replicate in treps[treps['assembly'] == ori_assembly(wildcards.assembly)].index:
            for function in [func for func in quality_control if
                             func.__name__ not in ['get_peak_calling_qc', 'get_trimming_qc']]:
                qc['files'].update(function(replicate))
            # scatac seq only on treps, not on single samples
            # and fastp also on treps
            if config.get("trimmer") == "fastp" or (get_workflow() == "scatac_seq" and get_trimming_qc in quality_control):
                qc['files'].update(get_trimming_qc(replicate))

    # qc on combined biological replicates/samples
    if get_peak_calling_qc in quality_control:
        for trep in treps[treps['assembly'] == ori_assembly(wildcards.assembly)].index:
            qc['files'].update(get_peak_calling_qc(trep))

    if get_rna_qc in quality_control:
        # Skip if desired, or if no BAMs are made (no trackhub + Salmon)
        if not (config.get("ignore_strandedness") or (config.get('quantifier', '') == 'salmon' and config.get('create_trackhub') == False)):
            for trep in treps[treps['assembly'] == ori_assembly(wildcards.assembly)].index:
                qc['files'].update(get_rna_qc(trep))

        # add dupRadar plots if BAMs are made
        if "REMOVE_DUPLICATES=true" not in config.get("markduplicates","") and not (config.get('quantifier', '') == 'salmon' and config.get('create_trackhub') == False):
            qc['files'].update(expand("{qc_dir}/dupRadar/{{assembly}}-dupRadar_mqc.png",**config))

    # DESeq2 sample distance/correlation cluster heatmaps
    if (get_peak_calling_qc in quality_control or get_rna_qc in quality_control) and len(treps.index) > 2:
        plots = ["sample_distance_clustering", "pearson_correlation_clustering", "spearman_correlation_clustering"]
        files = expand("{qc_dir}/plotCorrelation/{{assembly}}-DESeq2_{plots}_mqc.png", plots=plots, **config)
        # only perform clustering if there are 2 or more groups in the assembly
        for assembly in treps.assembly:
            if len(treps[treps.assembly == assembly].index) < 2:
                files = [f for f in files if f"clustering/{assembly}-" not in f]
        qc['files'].update(files)

    return qc


rule combine_qc_files:
    input:
        unpack(get_qc_files)
    output:
        expand("{qc_dir}/multiqc_{{assembly}}.tmp.files", **config),
    run:
        with open(output[0], mode="w") as out:
            out.write('\n'.join(input.files))


def get_qc_schemas(wildcards):
    qc = dict()
    qc['header'] = expand('{qc_dir}/header_info.yaml', **config)[0]
    qc['sample_names'] = expand('{qc_dir}/sample_names_{{assembly}}.tsv', **config)[0]
    qc['schema'] = expand('{qc_dir}/schema.yaml', **config)[0]
    if config["trimmer"] == "trimgalore":
        qc['filter_buttons'] = expand('{qc_dir}/sample_filters_{{assembly}}.tsv', **config)[0]
    return qc


rule multiqc:
    """
    Aggregate all the quality control metrics for every sample into a single 
    multiqc report.
    
    The input can get very long (causing problems with the shell), so we have 
    to write the input to a file, and then we can use that file as input to 
    multiqc.
    """
    input:
        unpack(get_qc_schemas),
        tmp=expand('{qc_dir}/samplesconfig_mqc.html', **config),
        files=rules.combine_qc_files.output,
    output:
        expand("{qc_dir}/multiqc_{{assembly}}.html", **config),
        directory(expand("{qc_dir}/multiqc_{{assembly}}_data", **config))
    message: explain_rule("multiqc")
    params:
        dir = "{qc_dir}/".format(**config),
        fqext1 = '_' + config['fqext1'],
        filter_buttons=lambda wildcards, input: f"--sample-filters {input.filter_buttons}" if hasattr(input, "filter_buttons") else ""
    log:
        expand("{log_dir}/multiqc_{{assembly}}.log", **config)
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc $(< {input.files}) -o {params.dir} -n multiqc_{wildcards.assembly}.html \
        --config {input.schema}                                                    \
        --config {input.header}                                                    \
        --sample-names {input.sample_names}                                        \
        {params.filter_buttons}                                    \
        --cl_config "extra_fn_clean_exts: [                                        \
            {{'pattern': ^.*{wildcards.assembly}-, 'type': 'regex'}},              \
            {{'pattern': {params.fqext1},          'type': 'regex'}},              \
            {{'pattern': _allsizes,                'type': 'regex'}},              \
            ]" > {log} 2>&1
        """


def get_trimming_qc(sample):
    if config["trimmer"] == "trimgalore":
        if get_workflow() == "scatac_seq":
            # we (at least for now) do not want fastqc for each single cell before and after trimming (too much for MultiQC).
            # still something to think about to add later, since that might be a good quality check though.
            if sampledict[sample]['layout'] == 'SINGLE':
                return expand(f"{{qc_dir}}/fastqc/{sample}_trimmed_fastqc.zip", **config)
            else:
                return expand(f"{{qc_dir}}/fastqc/{sample}_{{fqext}}_trimmed_fastqc.zip", **config)
        elif get_workflow() == "scrna_seq":
            # single-cell RNA seq does weird things with barcodes in the fastq file
            # therefore we can not just always start trimming paired-end even though
            # the samples are paired-end (ish)
            read_id = get_bustools_rid(config.get("count"))
            if read_id == 0:
                return expand([f"{{qc_dir}}/fastqc/{sample}_{{fqext1}}_fastqc.zip",
                               f"{{qc_dir}}/fastqc/{sample}_{{fqext1}}_trimmed_fastqc.zip",
                               f"{{qc_dir}}/trimming/{sample}_{{fqext1}}.{{fqsuffix}}.gz_trimming_report.txt"],
                              **config)
            elif read_id == 1:
                return expand([f"{{qc_dir}}/fastqc/{sample}_{{fqext2}}_fastqc.zip",
                            f"{{qc_dir}}/fastqc/{sample}_{{fqext2}}_trimmed_fastqc.zip",
                            f"{{qc_dir}}/trimming/{sample}_{{fqext2}}.{{fqsuffix}}.gz_trimming_report.txt"],
                            **config)
            else:
                raise NotImplementedError
        else:
            if sampledict[sample]['layout'] == 'SINGLE':
                return expand([f"{{qc_dir}}/fastqc/{sample}_fastqc.zip",
                               f"{{qc_dir}}/fastqc/{sample}_trimmed_fastqc.zip",
                               f"{{qc_dir}}/trimming/{sample}.{{fqsuffix}}.gz_trimming_report.txt"],
                               **config)
            else:
                return expand([f"{{qc_dir}}/fastqc/{sample}_{{fqext}}_fastqc.zip",
                               f"{{qc_dir}}/fastqc/{sample}_{{fqext}}_trimmed_fastqc.zip",
                               f"{{qc_dir}}/trimming/{sample}_{{fqext}}.{{fqsuffix}}.gz_trimming_report.txt"],
                               **config)

    elif config["trimmer"] == "fastp":
        if get_workflow() == "scrna_seq": 
            read_id = get_bustools_rid(config.get("count"))
            if read_id == 0:
                return expand(f"{{qc_dir}}/trimming/{sample}_{{fqext1}}.fastp.json", **config)
            elif read_id == 1:
                return expand(f"{{qc_dir}}/trimming/{sample}_{{fqext2}}.fastp.json", **config)
            else:
                raise NotImplementedError
        # not sure how fastp should work with scatac here
        return expand(f"{{qc_dir}}/trimming/{sample}.fastp.json", **config)


def get_alignment_qc(sample):
    output = []

    # add samtools stats
    output.append(f"{{qc_dir}}/markdup/{{{{assembly}}}}-{sample}.samtools-coordinate.metrics.txt")
    output.append(f"{{qc_dir}}/samtools_stats/{{aligner}}/{{{{assembly}}}}-{sample}.samtools-coordinate.samtools_stats.txt")
    
    # if Salmon is used, the sieving does not effect the expression values, so adding it to the MultiQC is confusing
    if sieve_bam(config) and \
            not (get_workflow() == "rna_seq" and config.get('quantifier') == 'salmon'):
        output.append(f"{{qc_dir}}/samtools_stats/{os.path.basename(config['final_bam_dir'])}/{{{{assembly}}}}-{sample}.samtools-coordinate.samtools_stats.txt")

    # add insert size metrics
    if sampledict[sample]['layout'] == "PAIRED":
        # if we do any sieving, we use the
        if config.get("filter_on_size"):
            output.append(f"{{qc_dir}}/InsertSizeMetrics/{{{{assembly}}}}-{sample}_allsizes.tsv")
        else:
            output.append(f"{{qc_dir}}/InsertSizeMetrics/{{{{assembly}}}}-{sample}.tsv")

    # get the ratio mitochondrial dna
    output.append(f"{{result_dir}}/{config['aligner']}/{{{{assembly}}}}-{sample}.samtools-coordinate-unsieved.bam.mtnucratiomtnuc.json")

    if get_workflow() in ["alignment", "chip_seq", "atac_seq", "scatac_seq"] or \
            get_workflow() == "rna_seq" and (config.get('create_trackhub') or config.get('quantifier') != 'salmon'):
        output.append("{qc_dir}/plotFingerprint/{{assembly}}.tsv")
    if len(breps[breps["assembly"] == treps.loc[sample, "assembly"]].index) > 1 and config.get("deeptools_qc"):
        output.append("{qc_dir}/plotCorrelation/{{assembly}}-deepTools_pearson_correlation_clustering_mqc.png")
        output.append("{qc_dir}/plotCorrelation/{{assembly}}-deepTools_spearman_correlation_clustering_mqc.png")
        output.append("{qc_dir}/plotPCA/{{assembly}}.tsv")

    return expand(output, **config)


def get_rna_qc(sample):
    output = []

    # add infer experiment reports
    col = samples.technical_replicates if "technical_replicates" in samples else samples.index
    if "strandedness" not in samples or samples[col == sample].strandedness[0] == "nan":
        output += expand(f"{{qc_dir}}/strandedness/{{{{assembly}}}}-{sample}.strandedness.txt",**config)

    return output


def get_peak_calling_qc(sample):
    output = []

    # add frips score (featurecounts)
    output.extend(expand(f"{{qc_dir}}/{{peak_caller}}/{{{{assembly}}}}-{sample}_featureCounts.txt.summary", **config))

    # macs specific QC
    if "macs2" in config["peak_caller"]:
        output.extend(expand(f"{{result_dir}}/macs2/{{{{assembly}}}}-{sample}_peaks.xls", **config))

    # deeptools profile
    assembly = treps.loc[sample, "assembly"]
    narrowpeak_used = "narrowPeak" in [get_peak_ftype(pc) for pc in list(config["peak_caller"].keys())]
    # TODO: replace with genomepy checkpoint in the future
    if has_annotation(assembly):
        output.extend(expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config))  # added to be unzipped
        if config.get("deeptools_qc"):
            output.extend(expand("{qc_dir}/plotProfile_gene/{{assembly}}-{peak_caller}.tsv", **config))
        if narrowpeak_used:
            output.extend(expand("{qc_dir}/chipseeker/{{assembly}}-{peak_caller}_img1_mqc.png", **config))
            output.extend(expand("{qc_dir}/chipseeker/{{assembly}}-{peak_caller}_img2_mqc.png", **config))
    if config.get("deeptools_qc") and narrowpeak_used:
        output.extend(expand("{qc_dir}/plotHeatmap_peaks/N{heatmap_npeaks}-{{assembly}}-deepTools_{peak_caller}_heatmap_mqc.png", **config))

    return output
