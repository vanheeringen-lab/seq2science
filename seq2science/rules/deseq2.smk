from seq2science.util import parse_de_contrasts

def get_contrasts():
    """
    splits contrasts that contain multiple comparisons
    """
    if not config.get("contrasts"):
        return []

    new_contrasts = []
    for contrast in list(config["contrasts"]):
        parsed_contrast, batch = parse_de_contrasts(contrast)
        components = len(parsed_contrast)

        # check the contrast column
        column_name = parsed_contrast[0]
        if column_name not in samples:
            raise IndexError(f"the DESeq2 contrast '{contrast}'\nreferences a column not found in the samples file: '{column_name}'.")

        if components == 2 or (components == 3 and "all" in contrast[1:]):
            # all vs 1 comparison ("A" or "all vs A")
            reflvl = str(parsed_contrast[2]) if parsed_contrast[1] == "all" else str(parsed_contrast[1])
            lvls = [str(lvl) for lvl in samples[column_name].dropna().unique()]
            lvls.remove(reflvl)
        else:
            # 1 vs 1 comparison ("A vs B")
            reflvl = parsed_contrast[2]
            lvls = [parsed_contrast[1]]

        for lvl in lvls:
            new_contrast = f"{column_name}_{lvl}_{reflvl}"
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
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/deseq2.R"


rule blind_clustering:
    """
    Create a sample distance matrix plot per assembly
    """
    input:
        expand("{counts_dir}/{{assembly}}-counts.tsv", **config),
    output:
        expand("{qc_dir}/clustering/{{assembly}}-Sample_clustering_mqc.png", **config),
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
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/deseq2_clustering.R"
