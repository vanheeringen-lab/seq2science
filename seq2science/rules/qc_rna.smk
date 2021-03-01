rule dupRadar:
    """
    visualize fraction of artifactual reads to normal read duplication
    (due to natural over-sequencing of highly expressed genes).
    """
    input:
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam",**config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf",**config),
        required=_strandedness_report,
    output:
        expand("{qc_dir}/dupRadar/{{assembly}}-dupRadar_{{sample}}_mqc.png", **config),
    log:
        expand("{log_dir}/dupRadar/{{assembly}}-{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/dupRadar/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    message: explain_rule("dupradar")
    params:
        strandedness=lambda wildcards: strandedness_to_quant(wildcards,"featurecounts"),
        paired=lambda wildcards: sampledict[wildcards.sample]["layout"] == "PAIRED",
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
        mem_gb=4,
    threads: 4
    conda:
        "../envs/dupradar.yaml"
    script:
        f"{config['rule_dir']}/../scripts/dupradar.R"
