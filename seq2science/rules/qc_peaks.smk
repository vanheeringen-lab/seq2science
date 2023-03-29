"""
all rules/logic related to peak calling quality control should be here.
"""


# TODO: featurecounts was not set up to accept hmmratac's gappedPeaks,
#  I added this as it was sufficient for rule idr
def get_featureCounts_bam(wildcards):
    if wildcards.peak_caller in ["macs2", "hmmratac"]:
        return expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    return expand("{final_bam_dir}/{{assembly}}-{{sample}}.sambamba-queryname.bam", **config)


# TODO: featurecounts was not set up to accept hmmratac's gappedPeaks,
#  I added this as it was sufficient for rule idr
def get_featureCounts_peak(wildcards):
    peak_sample = BREP_FROM_TREPS[wildcards.sample]
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
        summary=expand("{qc_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_featureCounts.txt.summary", **config),
    message: EXPLAIN["featureCounts_qc"]
    log:
        expand("{log_dir}/featureCounts/{{assembly}}-{{sample}}-{{peak_caller}}.log", **config),
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


rule plotHeatmap_peak:
    """
    Plot the peak heatmap using deepTools.
    """
    input:
        rules.computeMatrix_peak.output,
    output:
        img=expand(
            "{qc_dir}/plotHeatmap_peaks/N{{nrpeaks}}-{{assembly}}-deepTools_{{peak_caller}}_heatmap_mqc.png", **config
        ),
    log:
        expand("{log_dir}/plotHeatmap_peaks/{{assembly}}-{{peak_caller}}_N{{nrpeaks}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/plotHeatmap_peaks/{{assembly}}-{{peak_caller}}_N{{nrpeaks}}.benchmark.txt", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    resources:
        deeptools_limit=lambda wildcards, threads: threads,
    params:
        params=config.get("deeptools_heatmap_options", ""),
        slop=config.get("heatmap_slop", 0),
    shell:
        """
        plotHeatmap --matrixFile {input} --outFileName {output.img} {params.params} --startLabel="-{params.slop}" --endLabel={params.slop} > {log} 2>&1
        """


def get_summits_bed(wildcards):
    return expand(
        [
            f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_summits.bed"
            for replicate in breps[breps["assembly"] == ORI_ASSEMBLIES[wildcards.assembly]].index
        ],
        **config,
    )


rule chipseeker:
    """
    Make chipseeker plots, with the percentage of peaks in e.g. promoters.
    """
    input:
        narrowpeaks=get_summits_bed,
        gtf=rules.get_genome_annotation.output.gtf,
    output:
        img1=expand("{qc_dir}/chipseeker/{{assembly}}-{{peak_caller}}_img1_mqc.png", **config),
        img2=expand("{qc_dir}/chipseeker/{{assembly}}-{{peak_caller}}_img2_mqc.png", **config),
    params:
        names=lambda wildcards, input: get_descriptive_names(wildcards, input.narrowpeaks),
    log:
        expand("{log_dir}/chipseeker/{{assembly}}-{{peak_caller}}.log", **config),
    conda:
        "../envs/chipseeker.yaml"
    message: EXPLAIN["chipseeker"]
    resources:
        R_scripts=1,  # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/chipseeker.R"


rule upset_plot_peaks:
    """
    make an upset plot which shows the distribution of peaks amongst the biological replicates.
    """
    input:
        expand("{counts_dir}/{{peak_caller}}/{{assembly}}_onehotpeaks.tsv", **config)
    output:
        expand("{qc_dir}/upset/{{assembly}}-{{peak_caller}}_upset_mqc.jpg", **config)
    conda:
        "../envs/upset.yaml"
    script:
        f"{config['rule_dir']}/../scripts/upset.py"


rule maelstrom_report_preparation:
    """
    This rule injects the symlinked motif logos into the html so they get rendered in the multiqc report.
    """
    input:
        rules.gimme_maelstrom.output
    output:
        expand("{qc_dir}/gimme/{{assembly}}-{{gimme_database}}-{{peak_caller}}_mqc.html", **config)
    log:
        expand("{log_dir}/maelstrom_report_preparation/{{assembly}}-{{gimme_database}}-{{peak_caller}}.log", **config),
    script:
        f"{config['rule_dir']}/../scripts/maelstrom_report.py"
