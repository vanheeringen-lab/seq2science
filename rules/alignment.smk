config['bwaindex_types'] = ['amb', 'ann', 'bwt', 'pac', 'sa']


rule bwa_index:
    input:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
    output:
        expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
    log:
        expand("{log_dir}/bwa_index/{{assembly}}.log", **config)
    params:
        prefix="{genome_dir}/{{assembly}}/index/bwa/{{assembly}}".format(**config),
        algorithm=config["bwa_algo"]
    conda:
        "../envs/alignment.yaml"
    shell:
        "bwa index -p {params.prefix} -a {params.algorithm} {input} > {log} 2>&1"


rule bwa_mem:
    input:
        reads=expand("{result_dir}/{trimmed_dir}/{{sample}}", **config),
        index=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
    output:
        pipe(expand("{result_dir}/{bwa_dir}/{{sample}}-{{assembly}}.bampipe", **config))
    log:
        expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}.log", **config)
    params:
        interleaved='true',
        index_dir=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}", **config),
    threads: 20
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        bwa mem -t {threads} {params.index_dir} {input.reads}/* 1> {output} 2> {log}
        """


rule sambamba_sort:
    input:
        rules.bwa_mem.output
    output:
        expand("{result_dir}/{bwa_dir}/{{sample}}-{{assembly}}.bam", **config)
    log:
        expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}-sambamba_sort.log", **config)
    params:
        "-n" if config['bam_sort_order'] == 'queryname' else ''
    threads: 3
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        sambamba view --nthreads {threads} -S -f bam  {input[0]} -o /dev/stdout  2> {log} |
        sambamba sort --nthreads {threads} {params.sort}   /dev/stdin -o {output[0]}  2> {log}
        """


rule mark_duplicates:
    input:
        expand("{result_dir}/{bwa_dir}/{{sample}}-{{assembly}}.bam", **config)
    output:
        bam=    expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.bam", **config),
        metrics=expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.metrics.txt", **config)
    log:
        expand("{log_dir}/mark_duplicates/{{sample}}-{{assembly}}.log", **config)
    params:
        config['markduplicates']
    conda:
        "../envs/alignment.yaml"
    shell:
        "picard MarkDuplicates {params} INPUT={input} "
        "OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1"


# rule samtools_stats:
#     input:
#         expand("{result_dir}/dedup/{{sample}}-{{condition}}-{{project}}-{{assembly}}.bam", **config)
#     output:
#         expand("{result_dir}/samtools_stats/{{sample}}-{{condition}}-{{project}}-{{assembly}}.txt", **config)
#     log:
#         expand("{log_dir}/samtools_stats/{{sample}}-{{condition}}-{{project}}-{{assembly}}.log", **config)
#     conda:
#         "../envs/alignment.yaml"
#     shell:
#         "samtools stats {input} 1> {output} 2> {log}"
