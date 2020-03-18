def get_peak_replicates(wildcards):
    ftype = 'narrowPeak' if wildcards.peak_caller in ['macs2', 'genrich'] else 'gappedPeak'
    return expand([f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_peaks.{ftype}"
        for replicate in treps[treps['assembly'] == wildcards.assembly].index], **config)

rule peak_union:
    input:
        get_peak_replicates
    output:
        expand("{result_dir}/{{peak_caller}}/{{assembly}}_peaks.bed", **config)
    log:
        expand("{log_dir}/peak_union/{{assembly}}-{{peak_caller}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/peak_union/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/bedtools.yaml"
    shell:
        "cat {input} | sort -k1,1 -k2,2n | mergeBed -i stdin > {output} 2> {log}"


def get_coverage_table_replicates(file_ext):
    def wrapped(wildcards):
        return expand([f"{{dedup_dir}}/{wildcards.assembly}-{replicate}.{wildcards.sorter}-{wildcards.sorting}.{file_ext}"
            for replicate in treps[treps['assembly'] == wildcards.assembly].index], **config)
    return wrapped


rule coverage_table:
    input:
        peaks=rules.peak_union.output,
        replicates=get_coverage_table_replicates('bam'),
        replicate_bai=get_coverage_table_replicates('bam.bai')
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/count_table_{{assembly}}.{{sorter}}-{{sorting}}.txt", **config)
    log:
        expand("{log_dir}/multicov/{{assembly}}-{{peak_caller}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/multicov/{{assembly}}-{{peak_caller}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    conda:
        "../envs/gimme.yaml"
    shell:
        "coverage_table -p {input.peaks} -d {input.replicates} > {output} 2> {log}"
