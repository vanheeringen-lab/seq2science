"""
all rules/logic related to rna-specific quality control should be here.
"""


rule dupRadar:
    """
    visualize fraction of artifactual reads to normal read duplication
    (due to natural over-sequencing of highly expressed genes).
    """
    input:
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        report=rules.infer_strandedness.output,
    output:
        expand("{qc_dir}/dupRadar/{{assembly}}-{{sample}}.png", **config),
    log:
        expand("{log_dir}/dupRadar/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/dupRadar/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    message:
        EXPLAIN.get("dupradar", "")
    params:
        strandedness=lambda wildcards, input: get_strandedness(input.report[0], fmt="fc"),
        paired=lambda wildcards: SAMPLEDICT[wildcards.sample]["layout"] == "PAIRED",
    resources:
        R_scripts=1,  # conda's R can have issues when starting multiple times
        mem_gb=1,
    threads: 4
    conda:
        "../envs/dupradar.yaml"
    script:
        f"{config['rule_dir']}/../scripts/dupradar.R"


def get_dupradar_images(wildcards):
    output = []
    for trep in treps[treps["assembly"] == ORI_ASSEMBLIES[wildcards.assembly]].index:
        output += expand(f"{{qc_dir}}/dupRadar/{{{{assembly}}}}-{trep}.png", **config)
    return output


rule dupRadar_combine:
    """
    Combine the individual images (so we can group them nicely in the MultiQC).
    """
    input:
        get_dupradar_images,
    output:
        expand("{qc_dir}/dupRadar/{{assembly}}-dupRadar_mqc.png", **config),
    log:
        expand("{log_dir}/dupRadar/combine_{{assembly}}.log", **config),
    params:
        # created by the first dupRadar rule
        good_example=expand("{qc_dir}/dupRadar/good_example.png", **config),
        bad_example=expand("{qc_dir}/dupRadar/bad_example.png", **config),
    conda:
        "../envs/imagemick.yaml"
    shell:
        """
        convert {params.good_example} {params.bad_example} {input} -append {output} 2> {log}
        """
