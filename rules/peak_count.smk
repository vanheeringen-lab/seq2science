def get_peak_replicates(wildcards):
    if 'condition' in samples and config.get('combine_replicates', False):
        # if we have conditions use those peaks
        return expand([f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{condition}_peaks.narrowPeak"
            for condition in set(samples[samples['assembly'] == wildcards.assembly]['condition'])], **config)
    # otherwise all the separate ones
    return expand([f"{{result_dir}}/{wildcards.peak_caller}/{wildcards.assembly}-{replicate}_peaks.narrowPeak"
        for replicate in samples[samples['assembly'] == wildcards.assembly].index], **config)


rule peak_union:
    input:
        get_peak_replicates
    output:
        expand("{result_dir}/{{peak_caller}}/{{assembly}}_peaks.bed", **config)
    log:
        expand("{log_dir}/peak_union/{{assembly}}-{{peak_caller}}.log", **config)
    benchmark:
        expand("{log_dir}/peak_union/{{assembly}}-{{peak_caller}}.benchmark.txt", **config)[0]
    conda:
        "../envs/bedtools.yaml"
    shell:
        "cat {input} | sort -k1,1 -k2,2n | mergeBed -i stdin > {output}"


def get_multicov_replicates(file_ext):
    def wrapped(wildcards):
        if 'condition' in samples and config.get('combine_replicates', '') == 'merge':
            # if replicates' fastqs are merged get the merged bam
            return expand([f"{{dedup_dir}}/{wildcards.assembly}-{replicate}.{wildcards.sorter}-{wildcards.sorting}.{file_ext}"
                for replicate in set(samples[samples['assembly'] == wildcards.assembly]['condition'])], **config)
        # otherwise all the separate ones
        return expand([f"{{dedup_dir}}/{wildcards.assembly}-{replicate}.{wildcards.sorter}-{wildcards.sorting}.{file_ext}"
            for replicate in samples[samples['assembly'] == wildcards.assembly].index], **config)
    return wrapped


rule multicov:
    input:
        peaks=rules.peak_union.output,
        replicates=get_multicov_replicates('bam'),
        replicate_bai=get_multicov_replicates('bam.bai')
    output:
        expand("{result_dir}/count_table/{{peak_caller}}/count_table_{{assembly}}.{{sorter}}-{{sorting}}.bed", **config)
    log:
        expand("{log_dir}/multicov/{{assembly}}-{{peak_caller}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{log_dir}/multicov/{{assembly}}-{{peak_caller}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    conda:
        "../envs/bedtools.yaml"
    params:
        files=lambda wildcards, input: "\\t".join([file.split('-')[-2].split('.')[0] for file in input.replicates])
    shell:
        """
        echo -e "QNAME\\tstart\\tend\\t{params.files}" > {output}
        bedtools multicov -bams {input.replicates} -bed {input.peaks} >> {output}
        """
