"""
all rules/logic specific to scRNA qc and post-processing should be specified here.
"""


def get_count_dir(wildcards):
    # Return quantifier specific output directory
    if config["quantifier"] == "kallistobus":
        return rules.kallistobus_count.output.dir[0]
    else:
        return rules.citeseqcount.output.dir[0]


rule export_sce_obj:
    """
    Read scRNA UMI counts into a SingleCellExperiment object, add colData and export to RData format.
    """
    input:
        counts=get_count_dir,
    output:
        dir=expand(
            "{result_dir}/scrna-preprocess/{quantifier}/export/{{assembly}}-{{sample}}/{file}",
            **{**config, **{"file": ["raw_sce_obj.RData", "raw_seu_obj.RData"]}}
        ),
    log:
        expand("{log_dir}/scrna-preprocess/{quantifier}/export/{{assembly}}-{{sample}}_raw_sce_export.log", **config),
    priority: 1
    conda:
        "../envs/sce.yaml"
    params:
        isvelo=lambda wildcards, input: True if "--workflow lamanno" in config.get("count", "") else False,
        iskite=lambda wildcards, input: True if "--workflow kite" in config.get("count", "") else False,
        iscite=lambda wildcards, input: True if config["quantifier"] == "citeseqcount" else False,
        sample=lambda wildcards, input: rep_to_descriptive(wildcards.sample),
        replicates=True if "technical_replicates" in samples else False,
        scripts_dir=f"{config['rule_dir']}/../scripts/deseq2",
        outdir=lambda wildcards, input, output: os.path.dirname(output[0]),
    message:
        explain_rule("sce")
    resources:
        R_scripts=1,  # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/singlecell/read_kb_counts.R"


rule sctk_qc:
    """
    Perform scRNA QC with singleCellTK, store output in SingleCellExperiment object and export to RData format.
    """
    input:
       rds_raw=rules.export_sce_obj.output.dir[0]
    output:
        dir=expand(
            "{result_dir}/scrna-preprocess/{quantifier}/sctk/{{assembly}}-{{sample}}/{file}",
            **{**config, **{"file": ["export/sctk_sce_obj.RData", "SCTK_CellQC_summary.csv", "SCTK_CellQC.html"]}}
        ),
    log:
        expand("{log_dir}/scrna-preprocess/{quantifier}/sctk/{{assembly}}-{{sample}}_sctk.log", **config),
    priority: 1
    conda:
        "../envs/sctk.yaml"
    threads: 4
    params:
        sample=lambda wildcards, input: rep_to_descriptive(wildcards.sample),
        outdir=lambda wildcards, input, output: os.path.dirname(output[1]),
        isvelo=lambda wildcards, input: True if "--workflow lamanno" in config.get("count", "") else False,
        replicates=True if "technical_replicates" in samples else False,
    message:
        explain_rule("sctk")
    resources:
        R_scripts=1,
        mem_gb=50,# conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/singlecell/sctk_qc.R"
