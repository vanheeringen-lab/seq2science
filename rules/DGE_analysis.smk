def get_contrasts():
    """
    splits contrasts that contain multiple comparisons
    """
    if not config.get('contrasts', False):
        return []

    # contrasts from config
    old_contrasts = list(config["contrasts"])

    new_contrasts = []
    for contrast in old_contrasts:
        original_contrast, contrast, batch = parse_DE_contrasts(contrast)

        l = len(contrast)
        if l == 1:
            # write out the full contrast
            lvls = samples[contrast[0]].dropna().unique().tolist()
            new_contrast = batch + '+' + contrast[0] + '_' + lvls[0] + '_' + lvls[1] if batch is not None else contrast[0] + '_' + lvls[0] + '_' + lvls[1]
            new_contrasts.append(new_contrast)
        elif l == 2 or (l == 3 and 'all' in contrast[1:]):
            # create a list of all contrasts designed by the 'groupA vs all' design
            reflvl = str(contrast[2]) if contrast[1] == 'all' else str(contrast[1])
            lvls = samples[contrast[0]].dropna().unique().tolist()
            lvls = list(map(lambda x: str(x), lvls))
            lvls.remove(reflvl)

            for lvl in lvls:
                new_contrast = batch + '+' + contrast[0] + '_' + lvl + '_' + reflvl if batch is not None else contrast[0] + '_' + lvl + '_'  + reflvl
                new_contrasts.append(new_contrast)
        else:
            # remove '~', for uniformity
            new_contrast = original_contrast.replace("~", "").replace(" ", "")
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

workarounds:
-the workflow can be rerun until all scripts have finished
-the max number of parallel R-scripts can be set to 1
"""
rule deseq2:
    """
    Differential gene expression analysis with DESeq2.
    """
    input:
        expand("{counts_dir}/{{assembly}}-counts.tsv", **config)
    output:
        expand("{dge_dir}/{{assembly}}-{{contrast}}.diffexp.tsv", **config)
    conda:
        "../envs/deseq2.yaml"
    log:
        expand("{log_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.log", **config)
    benchmark:
        expand("{benchmark_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.benchmark.txt", **config)[0]
    threads: 4
    params:
        os.path.abspath(config["samples"])
    resources:
        R_scripts=1 # conda's R can have issues when starting multiple times
    script:
        "../scripts/deseq2.R"


rule blind_clustering:
    """
    Create a sample distance matrix plot per assembly
    """
    input:
         expand("{counts_dir}/{{assembly}}-counts.tsv", **config)
    output:
         expand("{dge_dir}/{{assembly}}-clustering.svg", **config)
    conda:
        "../envs/deseq2.yaml"
    log:
        expand("{log_dir}/deseq2/{{assembly}}-clustering.log", **config)
    benchmark:
        expand("{benchmark_dir}/deseq2/{{assembly}}-clustering.benchmark.txt", **config)[0]
    threads: 4
    params:
        os.path.abspath(config["samples"])
    resources:
        R_scripts=1 # conda's R can have issues when starting multiple times
    script:
        "../scripts/deseq2_clustering.R"
