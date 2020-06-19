ruleorder: macs2_callpeak > narrowpeak_summit

def count_table_output():
    """
    get all the count table outputs.
    """
    ftype = get_ftype(list(config['peak_caller'].keys())[0])
    if ftype != "narrowPeak":
        return []

    return expand(["{result_dir}/count_table/{peak_caller}/{assemblies}_{normalizes}.samtools-coordinate.tsv"],
                  **{**config,
                     **{'assemblies': set(samples['assembly']),
                        'peak_caller': config['peak_caller'].keys(),
                        'normalizes': ['raw', 'log1p_2', 'log1p_e', 'log1p_10', "quantilenorm_2", "quantilenorm_e",
                                       "quantilenorm_10"]}})


def get_peakfile_for_summit(wildcards):
    ftype = get_ftype(wildcards.peak_caller)
    if ftype != "narrowPeak":
        raise NotImplementedError("Narrowpeak to summit conversion is not supported for anything other than narrowpeak files. "
                                  "This means that we do not support peak counts for broadpeaks & gappedpeak format. Turn off: TODO")
    return expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_peaks.narrowPeak", **config)


rule narrowpeak_summit:
    """
    Convert a narrowpeak file to a "macs2 summits" file.
    """
    input:
        get_peakfile_for_summit
    output:
        expand("{result_dir}/{{peak_caller}}/{{assembly}}-{{sample}}_summits.bed", **config)
    log:
        expand("{log_dir}/bedtools_slop/{{sample}}-{{assembly}}-{{peak_caller}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/bedtools_slop/{{sample}}-{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{ print $1,$2+$10,$2+$10+1,$4,$9; }}' {input} > {output} 2> {log}
        """


def get_summitfiles(wildcards):
    return expand([f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_summits.bed"
        for replicate in treps[treps['assembly'] == wildcards.assembly].index], **config)


rule combine_peaks:
    """
    Uses gimmemotifs' combine_peaks to 
    """
    input:
        summitfiles=get_summitfiles,
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        temp(expand("{result_dir}/{{peak_caller}}/{{assembly}}_combinedsummits.bed", **config))
    log:
        expand("{log_dir}/bedtools_slop/{{assembly}}-{{peak_caller}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/bedtools_slop/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/gimme.yaml"
    shell:
        """
        combine_peaks -i {input.summitfiles} -g {input.sizes} > {output} 2> {log}
        """


rule bedtools_slop:
    """
    
    """
    input:
        bedfile=rules.combine_peaks.output,
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        temp(expand("{result_dir}/{{peak_caller}}/{{assembly}}_combinedpeaks.bed", **config))
    log:
        expand("{log_dir}/bedtools_slop/{{assembly}}-{{peak_caller}}.log", **config)
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
        return expand([f"{{final_bam_dir}}/{wildcards.assembly}-{replicate}.{wildcards.sorter}-{wildcards.sorting}.{file_ext}"
            for replicate in treps[treps['assembly'] == wildcards.assembly].index], **config)
    return wrapped


rule coverage_table:
    """
    
    """
    input:
        peaks=rules.bedtools_slop.output,
        replicates=get_coverage_table_replicates('bam'),
        replicate_bai=get_coverage_table_replicates('bam.bai')
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_raw.{{sorter}}-{{sorting}}.tsv", **config)
    log:
        expand("{log_dir}/multicov/{{assembly}}-{{peak_caller}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/multicov/{{assembly}}-{{peak_caller}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    conda:
        "../envs/gimme.yaml"
    shell:
        """
        echo "# The number of reads under each peak" > {output} 
        coverage_table -p {input.peaks} -d {input.replicates} 2> {log} | grep -vE "^#" 2>> {log} |  
        awk 'BEGIN {{ FS = "@" }} NR==1{{gsub("{wildcards.assembly}-|.samtools-coordinate","",$0)}}; \
        {{print $0}}' >> {output}
        """


rule log_normalization:
    """

    """
    input:
        rules.coverage_table.output
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_log1p_{{base}}.{{sorter}}-{{sorting}}.tsv", **config)
    run:
        import pandas as pd
        import numpy as np

        # read the coverage table as input
        cov_table = pd.read_csv(str(input), comment='#', index_col=0, sep="\t")

        # get our base
        if wildcards.base == "e":
            base = np.e
        else:
            base = float(wildcards.base)

        # take log1p
        species_log1p = np.log(cov_table + 1) / np.log(base)

        # save the df
        species_log1p.to_csv(str(output), index=True, index_label=str, header=True, sep="\t")

        # prepend a comment with how we normalized
        open(str(output), "w").write(
            f"# The number of reads under each peak, log1p normalized with base {wildcards.base}\n" +
            open(str(output)).read()
        )


rule quantile_normalization:
    """

    """
    input:
        rules.log_normalization.output
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/{{assembly}}_quantilenorm_{{base}}.{{sorter}}-{{sorting}}.tsv", **config)
    run:
        import pandas as pd
        import numpy as np

        def quantileNormalize(df_input):

            df = df_input.copy()

            # ranking the values and sort median on axis=1
            dic = {}
            for col in df:
                dic.update({col : sorted(df[col])})
            sorted_df = pd.DataFrame(dic)
            rank = sorted_df.median(axis = 1).tolist()

            # sort
            for col in df:
                t = np.searchsorted(np.sort(df[col]), df[col])
                df[col] = [rank[i] for i in t]

            return df

        species = pd.read_csv(str(input), comment='#', index_col=0, sep="\t")
        species_qN = quantileNormalize(species)
        species_qN.to_csv(str(output), index = True, index_label = str, header=True, sep="\t")
