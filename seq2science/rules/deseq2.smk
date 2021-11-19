from seq2science.util import expand_contrasts


def get_contrasts():
    """
    list all diffexp.tsv files we expect
    """
    if not config.get('contrasts'):
        return []

    contrasts=expand_contrasts(samples, config)
    all_contrasts = list()

    for de_contrast in contrasts:
        # parse groups
        target, reference = de_contrast.split("_")[-2:]
        column = de_contrast[:de_contrast.find(f"_{target}_{reference}")]

        # parse column
        if "+" in column:
            column = column.split("+")[1]
            
        if column not in samples:
            backup_columns = {
                "technical_replicates": "_trep",  # is trep technically possible? You need multiple reps right?
                "biological_replicates": "_brep",
                "descriptive_name": "_dname"
            }
            column = backup_columns[column]

        for assembly in all_assemblies:
            groups = set(samples[samples.assembly == assembly][column].to_list())
            if target in groups and reference in groups:
                all_contrasts.append(f"{config['deseq2_dir']}/{assembly}-{de_contrast}.diffexp.tsv")
    return all_contrasts


"""
Current issue with R:
when starting, ldpaths is updated. During this process it is opened and deleted.
Simultaneous R scripts can incidentally crash because it is missing.

active issue: https://github.com/conda-forge/r-base-feedstock/issues/67
active PR: https://github.com/conda/conda/pull/8776
"""

def deseq_input(wildcards):
    if "rna" in get_workflow():
        return expand("{counts_dir}/{{assembly}}-counts.tsv", **config)
    elif "atac"  in get_workflow() or "chip" in get_workflow():
        # only uses a single peak caller ------------------------------------------v
        # TODO different peak callers can probably be supported with wildcard_constraint peak_caller (.*) <-- empty allowed
        return expand("{counts_dir}/{peak_caller}/{{assembly}}_raw.tsv", **config)[0],
    else:
        raise NotImplementedError


# TODO once fixed the resource R_scripts can be removed
rule deseq2:
    """
    Differential gene expression analysis with DESeq2.
    """
    input:
        deseq_input
    output:
        diffexp=expand("{deseq2_dir}/{{assembly}}-{{contrast}}.diffexp.tsv", **config),
        maplot=expand("{qc_dir}/deseq2/{{assembly}}-{{contrast}}.ma_plot.png", **config),
        volcanoplot=expand("{qc_dir}/deseq2/{{assembly}}-{{contrast}}.volcano_plot.png", **config),
        pcaplot=expand("{qc_dir}/deseq2/{{assembly}}-{{contrast}}.pca_plot_mqc.png", **config),
    conda:
        "../envs/deseq2.yaml"
    log:
        expand("{log_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.log", **config),
    message: explain_rule("deseq2")
    benchmark:
        expand("{benchmark_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.benchmark.txt", **config)[0]
    threads: 4
    params:
        samples=os.path.abspath(config["samples"]),
        replicates=True if "technical_replicates" in samples else False,
        salmon=config.get("quantifier","") == "salmon",
        scripts_dir=f"{config['rule_dir']}/../scripts/deseq2"
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
        mem_gb=4,
    script:
        f"{config['rule_dir']}/../scripts/deseq2/deseq2.R"


rule merge_volcano_ma:
    """
    Combine the volcano and maplot resulting of the deseq2 rule into a single figure.
    """
    input:
        maplot=rules.deseq2.output.maplot,
        volcanoplot=rules.deseq2.output.volcanoplot
    output:
        expand("{qc_dir}/deseq2/{{assembly}}-{{contrast}}.combined_ma_volcano_mqc.png", **config)
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
        deseq_input
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
        R_scripts=1, # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/deseq2/clustering.R"
