"""
all rules/logic
"""
import os


"""
Current issue with R:
when starting, ldpaths is updated. During this process it is opened and deleted.
Simultaneous R scripts can incidentally crash because it is missing.

active issue: https://github.com/conda-forge/r-base-feedstock/issues/67
active PR: https://github.com/conda/conda/pull/8776
"""


def deseq_input(wildcards):
    if "rna" in WORKFLOW:
        return expand("{counts_dir}/{{assembly}}-counts.tsv", **config)
    elif "atac" in WORKFLOW or "chip" in WORKFLOW:
        # only uses a single peak caller ------------------------------------------v
        # TODO different peak callers can probably be supported with wildcard_constraint peak_caller (.*) <-- empty allowed
        return (expand("{counts_dir}/{peak_caller}/{{assembly}}_raw.tsv", **config)[0],)
    else:
        logger.error(
            f"The workflow you are running ({WORKFLOW}) does not support deseq2. "
            "Please make an issue on github if this is unexpected behaviour!"
        )
        os._exit(1)  # noqa


# TODO once fixed the resource R_scripts can be removed
rule deseq2:
    """
    Differential gene expression analysis with DESeq2.
    """
    input:
        deseq_input,
    output:
        diffexp=expand("{deseq2_dir}/{{assembly}}-{{contrast}}.diffexp.tsv", **config),
        maplot=expand("{qc_dir}/deseq2/{{assembly}}-{{contrast}}.ma_plot.png", **config),
        volcanoplot=expand("{qc_dir}/deseq2/{{assembly}}-{{contrast}}.volcano_plot.png", **config),
        pcaplot=expand("{qc_dir}/deseq2/{{assembly}}-{{contrast}}.pca_plot_mqc.png", **config),
    conda:
        "../envs/deseq2.yaml"
    log:
        expand("{log_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.log", **config),
    message: EXPLAIN["deseq2"]
    benchmark:
        expand("{benchmark_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.benchmark.txt", **config)[0]
    threads: 4
    params:
        samples=os.path.abspath(config["samples"]),
        replicates=True if "technical_replicates" in samples else False,
        salmon=config.get("quantifier", "") == "salmon",
        scripts_dir=f"{config['rule_dir']}/../scripts/deseq2",
    resources:
        R_scripts=1,  # conda's R can have issues when starting multiple times
        mem_gb=4,
    script:
        f"{config['rule_dir']}/../scripts/deseq2/deseq2.R"


rule merge_volcano_ma:
    """
    Combine the volcano and maplot resulting of the deseq2 rule into a single figure.
    """
    input:
        maplot=rules.deseq2.output.maplot,
        volcanoplot=rules.deseq2.output.volcanoplot,
    output:
        expand("{qc_dir}/deseq2/{{assembly}}-{{contrast}}.combined_ma_volcano_mqc.png", **config),
    log:
        expand("{log_dir}/deseq2/combine_{{assembly}}-{{contrast}}_plots.log", **config),
    conda:
        "../envs/imagemick.yaml"
    shell:
        """
        convert {input.maplot} {input.volcanoplot} +append {output} 2> {log}
        """


rule blind_clustering:
    """
    Create a sample distance cluster heatmap, 
    and various correlation cluster heatmaps.

    A (blind) variance stabilizing transformation is applied to the counts, 
    to reduce the impact of genes with low expression levels. 
    """
    input:
        deseq_input,
    output:
        expand("{qc_dir}/plotCorrelation/{{assembly}}-DESeq2_sample_distance_clustering_mqc.png", **config),
        expand("{qc_dir}/plotCorrelation/{{assembly}}-DESeq2_pearson_correlation_clustering_mqc.png", **config),
        expand("{qc_dir}/plotCorrelation/{{assembly}}-DESeq2_spearman_correlation_clustering_mqc.png", **config),
    log:
        expand("{log_dir}/plotCorrelation/{{assembly}}-DESeq2_clustering.log", **config),
    conda:
        "../envs/deseq2.yaml"
    threads: 4
    params:
        samples=os.path.abspath(config["samples"]),
        replicates=True if "technical_replicates" in samples else False,
        scripts_dir=f"{config['rule_dir']}/../scripts/deseq2",
    resources:
        R_scripts=1,  # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/deseq2/clustering.R"
