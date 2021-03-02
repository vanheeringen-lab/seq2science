ruleorder: macs2_callpeak > narrowpeak_summit

def count_table_output():
    """
    get all the count table outputs.
    """
    ftype = get_ftype(list(config["peak_caller"].keys())[0])
    if ftype != "narrowPeak":
        return []

    return expand(
        ["{counts_dir}/{peak_caller}/{assemblies}_{normalization}.tsv",
         "{counts_dir}/{peak_caller}/{assemblies}_onehotpeaks.tsv"],
        **{
            **config,
            **{
                "assemblies": all_assemblies,
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
        expand("{log_dir}/narrowpeak_summit/{{sample}}-{{assembly}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/narrowpeak_summit/{{sample}}-{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{ print $1,$2+$10,$2+$10+1,$4,$9; }}' {input} > {output} 2> {log}
        """


def get_summitfiles(wildcards):
    return expand(
        [
            f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_summits.bed"
            for replicate in breps[breps["assembly"] == ori_assembly(wildcards.assembly)].index
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
        expand("{log_dir}/combine_peaks/{{assembly}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/combine_peaks/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/gimme.yaml"
    params:
        windowsize=2 * config["peak_windowsize"],
        reps=lambda wildcards, input: input  # help resolve changes in input files
    message: explain_rule("combine_peaks")
    shell:
        """
        combine_peaks -i {input.summitfiles} -g {input.sizes} \
        --window {params.windowsize} > {output} 2> {log}
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
    params:
        slop=config["slop"],
        reps=lambda wildcards, input: input  # help resolve changes in input files
    message: explain_rule("bed_slop")
    shell:
        """
        bedtools slop -i {input.bedfile} -g {input.sizes} -b {params.slop} | uniq > {output} 2> {log}
        """


def get_coverage_table_replicates(file_ext):
    def wrapped(wildcards):
        if wildcards.peak_caller == "macs2":
            return expand(
                [
                    f"{{final_bam_dir}}/{wildcards.assembly}-{replicate}.samtools-coordinate.{file_ext}"
                    for replicate in treps[treps["assembly"] == ori_assembly(wildcards.assembly)].index
                ],
                **config,
            )
        elif wildcards.peak_caller == "genrich":
            return expand(
                [
                    f"{{final_bam_dir}}/{wildcards.assembly}-{replicate}.sambamba-queryname.{file_ext}"
                    for replicate in treps[treps["assembly"] == ori_assembly(wildcards.assembly)].index
                ],
                **config,
            )
        else:
            raise NotImplementedError

    return wrapped


rule coverage_table:
    """
    Use gimmemotif's coverage_table to generate a count table with the nr of reads
    under each peak per sample.
    """
    input:
        peaks=rules.bedtools_slop.output,
        replicates=get_coverage_table_replicates("bam"),
        replicate_bai=get_coverage_table_replicates("bam.bai"),
    output:
        expand("{counts_dir}/{{peak_caller}}/{{assembly}}_raw.tsv", **config),
    log:
        expand("{log_dir}/coverage_table/{{assembly}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/coverage_table/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/gimme.yaml"
    resources:
        mem_gb=3,
    message: explain_rule("coverage_table")
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
        expand("{counts_dir}/{{peak_caller}}/{{assembly}}_quantilenorm.tsv", **config),
    log:
        expand("{log_dir}/quantile_normalization/{{assembly}}-{{peak_caller}}-quantilenorm.log", **config),
    benchmark:
        expand("{benchmark_dir}/quantile_normalization/{{assembly}}-{{peak_caller}}-quantilenorm.benchmark.txt", **config)[0]
    conda:
        "../envs/qnorm.yaml"
    threads: 4
    script:
        f"{config['rule_dir']}/../scripts/qnorm.py"


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
        expand("{counts_dir}/{{peak_caller}}/{{assembly}}_{{normalisation,(TMM|RLE|upperquartile)}}.tsv", **config),
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
        expand("{counts_dir}/{{peak_caller}}/{{assembly}}_{{normalisation}}.tsv", **config),
    output:
        expand("{counts_dir}/{{peak_caller}}/{{assembly}}_log{{base}}_{{normalisation}}.tsv", **config),
    run:
        import pandas as pd
        import numpy as np

        # read the coverage table as input
        cov_table = pd.read_csv(str(input), comment="#", index_col=0, sep="\t")

        # get our base
        base = float(wildcards.base)

        if cov_table.isin(["NA"]).any(axis=None):
            # return all NA if NA exists
            species_log = cov_table.copy()
            species_log[:] = "NA"
        else:
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
        expand("{counts_dir}/{{peak_caller}}/{{assembly}}_log{{base}}_{{normalisation}}.tsv", **config),
    output:
        expand("{counts_dir}/{{peak_caller}}/{{assembly}}_meancenter_log{{base}}_{{normalisation}}.tsv", **config),
    run:
        import pandas as pd

        # read the coverage table as input
        cov_table = pd.read_csv(str(input), comment="#", index_col=0, sep="\t")

        # take mean center
        if cov_table.isin(["NA"]).any(axis=None):
            # return all NA if NA exists
            cov_mc = cov_table.copy()
            cov_mc[:] = "NA"
        else:
            cov_mc = cov_table.subtract(cov_table.mean(axis=1), axis=0)

        # prepend a comment with how we normalized
        open(str(output), "w").write(
            f"# The number of reads under each peak, mean centered after log1p {wildcards.base} and {wildcards.normalisation} normalization\n"
            + cov_mc.to_csv(index=True, header=True, sep="\t")
        )


def get_all_narrowpeaks(wildcards):
    return [f"{config['result_dir']}/{{peak_caller}}/{{assembly}}-{replicate}_peaks.narrowPeak"
        for replicate in breps[breps['assembly'] == ori_assembly(wildcards.assembly)].index]


rule onehot_peaks:
    """
    Get onehot encodings of which peaks are found in which samples
    """
    input:
        narrowpeaks=get_all_narrowpeaks,
        combinedpeaks=rules.combine_peaks.output
    output:
        real=expand("{counts_dir}/{{peak_caller}}/{{assembly}}_onehotpeaks.tsv", **config),
        tmp=temp(expand("{counts_dir}/{{peak_caller}}/{{assembly}}_onehotpeaks.tsv.tmp", **config))
    conda:
        "../envs/bedtools.yaml"
    log:
        expand("{log_dir}/onehot_peaks/{{assembly}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/onehot_peaks/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    params:
        reps=lambda wildcards, input: input  # help resolve changes in input files
    shell:
        """
        awk '{{print $1":"$2"-"$3}}' {input.combinedpeaks} > {output.real} 2> {log}

        for brep in {input.narrowpeaks}
        do
            bedtools sort -i $brep 2> {log} |
            bedtools intersect -a {input.combinedpeaks} -b stdin -c 2> {log} |
            awk '{{print $4}}' 2> {log} |
            paste {output.real} - > {output.tmp}
            mv {output.tmp} {output.real}
        done
        
        echo -e "# onehot encoding of which condition contains which peaks\\n$(cat {output.real})" > {output.tmp} && cp {output.tmp} {output.real}
        """
