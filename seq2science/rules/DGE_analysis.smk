def get_contrasts():
    """
    splits contrasts that contain multiple comparisons
    """
    if "contrasts" not in config:
        return []

    new_contrasts = []
    for contrast in list(config["contrasts"]):
        parsed_contrast, batch = parse_de_contrasts(contrast)
        column_name = parsed_contrast[0]
        components = len(parsed_contrast)
        if components == 2 or (components == 3 and "all" in contrast[1:]):
            # create a list of all contrasts designed by the 'groupA vs all' design
            reflvl = str(parsed_contrast[2]) if parsed_contrast[1] == "all" else str(parsed_contrast[1])
            lvls = [str(lvl) for lvl in samples[column_name].dropna().unique().tolist()]
            lvls.remove(reflvl)
        else:
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


# TODO once fixed the resource R_scripts can be removed
rule deseq2:
    """
    Differential gene expression analysis with DESeq2.
    """
    input:
        expand("{counts_dir}/{{assembly}}-counts.tsv", **config),
    output:
        expand("{dge_dir}/{{assembly}}-{{contrast}}.diffexp.tsv", **config),
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
        replicates=True if "replicate" in samples else False,
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
        expand("{dge_dir}/{{assembly}}-clustering.svg", **config),
    conda:
        "../envs/deseq2.yaml"
    log:
        expand("{log_dir}/deseq2/{{assembly}}-clustering.log", **config),
    benchmark:
        expand("{benchmark_dir}/deseq2/{{assembly}}-clustering.benchmark.txt", **config)[0]
    threads: 4
    params:
        samples=os.path.abspath(config["samples"]),
        replicates=True if "replicate" in samples else False,
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/deseq2_clustering.R"
