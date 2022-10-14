"""
all rules/logic specific to scRNA qc and post-processing should be specified here.
"""
import os


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
            "{result_dir}/scrna-preprocess/{quantifier}/raw/{{assembly}}-{{sample}}/{file}",
            **{**config, **{"file": ["export/R/raw_SCE.RDS", "SCE_raw_summary.csv"]}}
        ),
    log:
        expand("{log_dir}/scrna-preprocess/{quantifier}/raw/{{assembly}}-{{sample}}_raw_sce.log", **config),
    priority: 1
    conda:
        "../envs/sce.yaml"
    params:
        isvelo=lambda wildcards, input: True if "--workflow lamanno" in config.get("count", "") else False,
        iskite=lambda wildcards, input: True if "--workflow kite" in config.get("count", "") else False,
        iscite=lambda wildcards, input: True if config["quantifier"] == "citeseqcount" else False,
        sample=lambda wildcards, input: rep_to_descriptive(wildcards.sample),
        replicates=True if "technical_replicates" in samples else False,
        scripts_dir=f"{config['rule_dir']}/../scripts",
        outdir=lambda wildcards, input, output: os.path.dirname(output[1]),
    message: EXPLAIN["sce"]
    resources:
        R_scripts=1,  # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/singlecell/read_kb_counts.R"


rule get_mt_genes:
    """
    Extract MT genes from gtf annotation for QC
    """
    input:
        gtf=rules.get_genome_annotation.output.gtf,
    output:
        dir=expand("{result_dir}/scrna-preprocess/{quantifier}/{{assembly}}-mt.txt", **config),
    log:
        expand("{log_dir}/scrna-preprocess/{quantifier}/{{assembly}}-mt.log", **config),
    script:
        f"{config['rule_dir']}/../scripts/genomepy/get_genome_mt.py"
    

rule sctk_qc:
    """
    Perform scRNA QC with singleCellTK, store output in SingleCellExperiment object and export to RData format.
    """
    input:
       rds_raw=rules.export_sce_obj.output.dir[0]
    output:
        dir=expand(
            "{result_dir}/scrna-preprocess/{quantifier}/sctk/{{assembly}}-{{sample}}/{file}",
            **{**config, **{"file": ["export/R/sctk_SCE.RDS", "SCE_sctk_summary.csv"]}}
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
        scripts_dir=f"{config['rule_dir']}/../scripts",
        replicates=True if "technical_replicates" in samples else False,
    message: EXPLAIN["sctk"]
    resources:
        R_scripts=1,
        mem_gb=50,# conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/singlecell/sctk_qc.R"
