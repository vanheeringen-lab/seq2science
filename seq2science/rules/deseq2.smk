from seq2science.util import parse_contrast

def get_contrasts():
    """
    splits contrasts that contain multiple comparisons
    """
    if not config.get("contrasts"):
        return []

    new_contrasts = []
    for contrast in list(config["contrasts"]):
        batch, column, target, reference = parse_contrast(contrast, samples, check=False)

        if target == "all":
            # all vs 1 comparison ("all vs A")
            targets = set(samples[column].dropna().astype(str))
            targets.remove(reference)
        else:
            # 1 vs 1 comparison ("A vs B")
            targets = [target]

        for target in targets:
            new_contrast = f"{column}_{target}_{reference}"
            if batch:
                new_contrast = f"{batch}+{new_contrast}"
            new_contrasts.append(new_contrast)

    # get unique elements
    new_contrasts = list(set(new_contrasts))
    return new_contrasts


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
        expand("{deseq2_dir}/{{assembly}}-{{contrast}}.diffexp.tsv", **config),
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
        replicates=True if "technical_replicate" in samples else False,
        salmon=config.get("quantifier","") == "salmon",
        scripts_dir=f"{config['rule_dir']}/../scripts/deseq2"
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
        mem_gb=4,
    script:
        f"{config['rule_dir']}/../scripts/deseq2/deseq2.R"


rule blind_clustering:
    """
    Create a sample distance cluster heatmap, 
    and various correlation cluster heatmaps.
    
    A (blind) variance stabilizing transformation is applied to the counts, 
    to reduce the impact of genes with low expression levels. 
    """
    input:
        deseq_input  # expand("{counts_dir}/{{assembly}}-counts.tsv", **config),
    output:
        expand("{qc_dir}/clustering/{{assembly}}-Sample_distance_clustering_mqc.png", **config),
        expand("{qc_dir}/clustering/{{assembly}}-Pearson_correlation_clustering_mqc.png", **config),
        expand("{qc_dir}/clustering/{{assembly}}-Spearman_correlation_clustering_mqc.png", **config),
    conda:
        "../envs/deseq2.yaml"
    log:
        expand("{log_dir}/deseq2/{{assembly}}-clustering.log", **config),
    benchmark:
        expand("{benchmark_dir}/deseq2/{{assembly}}-clustering.benchmark.txt", **config)[0]
    threads: 4
    params:
        samples=os.path.abspath(config["samples"]),
        replicates=True if "technical_replicate" in samples else False,
        scripts_dir=f"{config['rule_dir']}/../scripts/deseq2",
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/deseq2/clustering.R"
