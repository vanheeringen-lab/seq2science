def get_contrasts():
    """
    splits contrasts that contain multiple comparisons
    """
    # contrasts from config
    old_contrasts = list(config["contrasts"])

    new_contrasts = []
    for contrast in old_contrasts:
        original_contrast = contrast

        # remove whitespaces (and '~'s if used)
        contrast = contrast.replace(" ", "")
        contrast = contrast.replace("~", "")

        # split and store batch effect
        batch = None
        if '+' in contrast:
            batch = contrast.split('+')[0]
            contrast = contrast.split('+')[1]

        # parse contrast
        contrast = contrast.split('_')

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
                new_contrast = batch + '+' + contrast[0] + '_' + reflvl + '_' + lvl if batch is not None else contrast[0] + '_' + reflvl + '_' + lvl
                new_contrasts.append(new_contrast)
        else:
            # remove '~', for uniformity
            new_contrast = original_contrast.replace("~", "").replace(" ", "")
            new_contrasts.append(new_contrast)

    # get unique elements
    new_contrasts = list(set(new_contrasts))

    return new_contrasts


rule deseq2:
    """
    Differential gene expression analysis with DESeq2.
    """
    input:
        expand("{result_dir}/gene_counts/{{assembly}}-counts.tsv", **config)
    output:
        table=expand("{result_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.tsv", **config),
        ma_plot=expand("{result_dir}/deseq2/{{assembly}}-{{contrast}}.ma_plot.svg", **config),
        pca_plot=expand("{result_dir}/deseq2/{{assembly}}-{{contrast}}.pca_plot.svg", **config)
    conda:
        "../envs/deseq2.yaml"
    log:
        expand("{log_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.log", **config)
    benchmark:
        expand("{benchmark_dir}/deseq2/{{assembly}}-{{contrast}}.diffexp.benchmark.txt", **config)[0]
    threads: 4
    params:
        os.path.abspath(config["samples"])
    script:
        "../scripts/deseq2.R"

# load DE analysis contrasts
contrasts = ['']
if config.get('contrasts', False):
    contrasts = get_contrasts()

rule blind_clustering:
    """
    Create a sample distance matrix plot per assembly
    """
    input:
         expand("{result_dir}/gene_counts/{{assembly}}-counts.tsv", **config)
    output:
         expand("{result_dir}/deseq2/{{assembly}}-clustering.svg", **config)
    conda:
        "../envs/deseq2.yaml"
    log:
        expand("{log_dir}/deseq2/{{assembly}}-clustering.log", **config)
    benchmark:
        expand("{benchmark_dir}/deseq2/{{assembly}}-clustering.benchmark.txt", **config)[0]
    threads: 4
    params:
        samples=os.path.abspath(config["samples"]),
        some_contrast=contrasts[0]
    script:
        "../scripts/deseq2_clustering.R"
