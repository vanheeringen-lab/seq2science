ruleorder: macs2_callpeak > narrowpeak_summit


def count_table_output():
    """
    get all the count table outputs.
    """
    ftype = get_ftype(list(config["peak_caller"].keys())[0])
    if ftype != "narrowPeak":
        return []

    return expand(
        ["{result_dir}/count_table/{peak_caller}/{assemblies}_{normalization}.tsv"],
        **{
            **config,
            **{
                "assemblies": set(samples["assembly"]),
                "peak_caller": config["peak_caller"].keys(),
                "normalization": [
                    "raw",
                    f"meancenter_log{config['logbase']}_quantilenorm",
                    f"meancenter_log{config['logbase']}_TMM",
                    f"meancenter_log{config['logbase']}_RLE",
                    f"meancenter_log{config['logbase']}_upperquartile",
                ],
                "mc": ["", "meancenter_"],
            },
        },
    )


def get_peakfile_for_summit(wildcards):
    ftype = get_ftype(wildcards.peak_caller)
    if ftype != "narrowPeak":
        raise NotImplementedError(
            "Narrowpeak to summit conversion is not supported for anything other than narrowpeak files. "
            "This means that we do not support peak counts for broadpeaks & gappedpeak format. Turn off: TODO"
        )
    return expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_peaks.narrowPeak", **config)


rule narrowpeak_summit:
    """
    Convert a narrowpeak file to a "macs2 summits" file.
    """
    input:
        get_peakfile_for_summit,
    output:
        expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_summits.bed", **config),
    log:
        expand("{log_dir}/bedtools_slop/{{sample}}-{{assembly}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bedtools_slop/{{sample}}-{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{ print $1,$2+$10,$2+$10+1,$4,$9; }}' {input} > {output} 2> {log}
        """


def get_summitfiles(wildcards):
    return expand(
        [
            f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_summits.bed"
            for replicate in treps[treps["assembly"] == wildcards.assembly].index
        ],
        **config,
    )


rule combine_peaks:
    """
    Uses gimmemotifs' combine_peaks to "combine" peaks. This finds all peaks close
    together and takes the most significant one as the true peak.
    """
    input:
        summitfiles=get_summitfiles,
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        temp(expand("{result_dir}/{{peak_caller}}/{{assembly}}_combinedsummits.bed", **config)),
    log:
        expand("{log_dir}/bedtools_slop/{{assembly}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bedtools_slop/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/gimme.yaml"
    params:
        2 * config["peak_windowsize"],
    shell:
        """
        combine_peaks -i {input.summitfiles} -g {input.sizes} \
        --window {params} > {output} 2> {log}
        """


rule bedtools_slop:
    """
    After combine_peaks we end up with just a bed file of summits. We extend all peaks
    to a total width of 200, for a fair comparison between peaks.
    """
    input:
        bedfile=rules.combine_peaks.output,
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
    output:
        temp(expand("{result_dir}/{{peak_caller}}/{{assembly}}_combinedpeaks.bed", **config)),
    log:
        expand("{log_dir}/bedtools_slop/{{assembly}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bedtools_slop/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools slop -i {input.bedfile} -g {input.sizes} -b {config[slop]} > {output} 2> {log}
        """


def get_coverage_table_replicates(file_ext):
    def wrapped(wildcards):
        if wildcards.peak_caller == "macs2":
            return expand(
                [
                    f"{{final_bam_dir}}/{wildcards.assembly}-{replicate}.samtools-coordinate.{file_ext}"
                    for replicate in treps[treps["assembly"] == wildcards.assembly].index
                ],
                **config,
            )
        elif wildcards.peak_caller == "genrich":
            return expand(
                [
                    f"{{final_bam_dir}}/{wildcards.assembly}-{replicate}.sambamba-queryname.{file_ext}"
                    for replicate in treps[treps["assembly"] == wildcards.assembly].index
                ],
                **config,
            )
        else:
            raise NotImplementedError

    return wrapped


rule coverage_table:
    """
    Use gimmemotif's coverage_table to generate a cpunt table with the nr of reads
    under each peak per sample.
    """
    input:
        peaks=rules.bedtools_slop.output,
        replicates=get_coverage_table_replicates("bam"),
        replicate_bai=get_coverage_table_replicates("bam.bai"),
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_raw.tsv", **config),
    log:
        expand("{log_dir}/multicov/{{assembly}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/multicov/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/gimme.yaml"
    shell:
        """
        echo "# The number of reads under each peak" > {output} 
        coverage_table -p {input.peaks} -d {input.replicates} 2> {log} | grep -vE "^#" 2>> {log} |  
        awk 'BEGIN {{ FS = "@" }} NR==1{{gsub("{wildcards.assembly}-|.samtools-coordinate","",$0)}}; \
        {{print $0}}' >> {output}
        """


rule quantile_normalization:
    """
    Quantile normalization is a type of normalization that makes the distribution 
    between samples identical. This means that the actual count distribution 
    within a sample changes. After quantile normalization, samples are CPM 
    normalized (which has no effect between and within samples), for a fair
    comparison between normalisation techniques.

    See wikipedia:
    https://en.wikipedia.org/wiki/Quantile_normalization
    """
    input:
        rules.coverage_table.output,
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_quantilenorm.tsv", **config),
    run:
        import pandas as pd


        def quantileNormalize_cpm(df):
            """
        solution adapted from:
        https://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe/41078786#41078786

        Changes:
        - takes the mean of ties within a sample instead of ignoring it (as per wiki example)
        - normalize on counts per million (cpm)

        Note: naive implementation, might get slow with large dataframes
        """
            # get all the ranks, where ties share a spot
            rank_ties = df.rank(method="min").astype(int)

            # calculate the median per rank
            rank_means = (df * 1_000_000 / df.sum(axis=0)).stack().groupby(df.rank(method="first").stack().astype(int)).mean()

            # quantile normalize and ignore ties
            qn_df = df.rank(method="min").stack().astype(int).map(rank_means).unstack()

            # now fix our ties by taking their mean
            for column in qn_df.columns:
                for idx, count in zip(*np.unique(qn_df[column].rank(method="min").astype(int), return_counts=True)):
                    if count > 1:
                        wrong_idxs = np.where(rank_ties[column] == idx)[0]
                        qn_df[column].iloc[wrong_idxs] = np.mean(rank_means[idx - 1 : idx - 1 + count])

            return qn_df


        df = pd.read_csv(str(input), comment="#", index_col=0, sep="\t")
        df_qn = quantileNormalize_cpm(df)
        open(str(output), "w").write(
            "# The number of reads under each peak, cpm quantile normalized\n"
            + df_qn.to_csv(index_label="loc", index=True, header=True, sep="\t")
        )



rule edgeR_normalization:
    """
    edgeR supports three different types of normalization: TMM, RLE, and upperquartile.

    TMM: is the weighted trimmed mean of M-values proposed by Robinson and Oshlack (2010).
    RLE: is the scaling factor method proposed by Anders and Huber (2010). DEseq2 normalisation
         is based on this.
    upperquartile: is the upper-quartile normalization method of Bullard et al (2010).

    In addition to 
    """
    input:
        rules.coverage_table.output,
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_{{normalisation,(TMM|RLE|upperquartile)}}.tsv", **config),
    log:
        expand("{log_dir}/edgeR_normalization/{{assembly}}-{{peak_caller}}-{{normalisation}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/edgeR_normalization/{{assembly}}-{{peak_caller}}-{{normalisation}}.benchmark.txt", **config)[0]
    conda:
        "../envs/edger.yaml"
    resources:
        R_scripts=1, # conda's R can have issues when starting multiple times
    script:
        f"{config['rule_dir']}/../scripts/edger_norm.R"


rule log_normalization:
    """
    Log1p normalization of a count table.
    """
    input:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_{{normalisation}}.tsv", **config),
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_log{{base}}_{{normalisation}}.tsv", **config),
    run:
        import pandas as pd
        import numpy as np

        # read the coverage table as input
        cov_table = pd.read_csv(str(input), comment="#", index_col=0, sep="\t")

        # get our base
        base = float(wildcards.base)

        # take log
        species_log = np.log1p(cov_table) / np.log(base)

        # prepend a comment with how we normalized
        open(str(output), "w").write(
            f"# The number of reads under each peak, log1p normalized with base {wildcards.base} after "
            f"{wildcards.normalisation} normalisation\n" + species_log.to_csv(index=True, header=True, sep="\t")
        )



rule mean_center:
    """
    Mean centering of a count table.
    """
    input:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_log{{base}}_{{normalisation}}.tsv", **config),
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_meancenter_log{{base}}_{{normalisation}}.tsv", **config),
    run:
        import pandas as pd
        import numpy as np

        # read the coverage table as input
        cov_table = pd.read_csv(str(input), comment="#", index_col=0, sep="\t")

        # take mean center
        cov_mc = cov_table.subtract(cov_table.mean(axis=1), axis=0)

        # prepend a comment with how we normalized
        open(str(output), "w").write(
            f"# The number of reads under each peak, mean centered after log1p {wildcards.base} and {wildcards.normalisation} normalization\n"
            + cov_mc.to_csv(index=True, header=True, sep="\t")
        )
