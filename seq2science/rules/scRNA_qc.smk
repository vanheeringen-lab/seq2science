"""
all rules/logic specific to scRNA qc and post-processing should be specified here.
"""


def get_count_dir(wildcards):
    # Return quantifier specific output directory
    if config["quantifier"] == "kallistobus":
        return rules.kallistobus_count.output.dir[0]
    else:
        return rules.citeseqcount.output.dir[0]


rule export_seurat_obj:
    """
    Read scRNA count output into Seurat object, add meta-data and export to RData format.
    """
    input:
        counts=get_count_dir,
    output:
        rds=expand("{result_dir}/seurat/{quantifier}/{{assembly}}-{{sample}}_seu_obj.RData", **config),
    log:
        expand("{log_dir}/seurat/{{assembly}}-{{sample}}_seu_obj.log", **config),
    priority: 1
    conda:
        "../envs/seurat.yaml"
    params:
        isvelo=lambda wildcards, input: True if "--workflow lamanno" in config.get("count", "") else False,
        iskite=lambda wildcards, input: True if "--workflow kite" in config.get("count", "") else False,
        iscite=lambda wildcards, input: True if config["quantifier"] == "citeseqcount" else False,
        sample=lambda wildcards, input: rep_to_descriptive(wildcards.sample),
        replicates=True if "technical_replicates" in samples else False,
        scripts_dir=f"{config['rule_dir']}/../scripts/deseq2",
    message:
        explain_rule("seurat")
    resources:
        R_scripts=1,  # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/seurat/read_kb_counts.R"


rule sctk_qc:
    """
    Read scRNA count output into Seurat object, add meta-data and export to RData format.
    """
    input:
       rds_raw=rules.export_seurat_obj.output.rds
    output:
        dir=expand(
            "{result_dir}/sctk/{quantifier}/{{assembly}}-{{sample}}/{file}",
            **{**config, **{"file": ["seu_obj_sctk.RData", "SCTK_CellQC_summary.csv", "SCTK_CellQC.html"]}}
        ),
    log:
        expand("{log_dir}/sctk/{{assembly}}-{{sample}}_sctk.log", **config),
    priority: 1
    conda:
        "../envs/sctk.yaml"
    params:
        sample=lambda wildcards, input: rep_to_descriptive(wildcards.sample),
        outdir=lambda wildcards, input, output: os.path.dirname(output[0]),
        isvelo=lambda wildcards, input: True if "--workflow lamanno" in config.get("count", "") else False,
        replicates=True if "technical_replicates" in samples else False,
    message:
        explain_rule("sctk")
    resources:
        R_scripts=1,
        mem_gb=50,# conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/seurat/sctk_qc.R"


def get_merge_objs(wildcards):
    # Return quantifier specific output directory
    if config["run_sctk_qc"]:
        return expand(
            [
                f"{{result_dir}}/sctk/{{quantifier}}/{custom_assembly(treps.loc[trep, 'assembly'])}-{trep}/seu_obj_sctk.RData"
                for trep in treps.index
            ],
            **config,
        )
    else:
        return expand(
            [
                f"{{result_dir}}/seurat/{{quantifier}}/{custom_assembly(treps.loc[trep, 'assembly'])}-{trep}_seu_obj.RData"
                for trep in treps.index
            ],
            **config,
        )


rule merge_seurat_obj:
    """
    Gather and merge multiple Seurat objects into a combined object and export to RData format.
    """
    input:
        seu_objs=get_merge_objs,
    output:
        rds=f"{config['result_dir']}/seurat/{{quantifier}}/{{assembly}}_seu_merged.RData",
    log:
        expand("{log_dir}/seurat/{{quantifier}}/{{assembly}}_seu_merged.log", **config),
    priority: 1
    conda:
        "../envs/seurat.yaml"
    params:
        isvelo=lambda wildcards, input: True if "--workflow lamanno" in config.get("count", "") else False,
    message:
        explain_rule("seurat")
    resources:
        R_scripts=1,  # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/seurat/merge_seurat_objs.R"
