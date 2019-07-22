config['bwaindex_types'] = ['amb', 'ann', 'bwt', 'pac', 'sa']

# TODO replace wrapper with one-liner?
# TODO merge params
rule bwa_index:
    input:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
    output:
        expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
    log:
        "logs/bwa_index/{assembly}.log"
    params:
        prefix="{genome_dir}/{{assembly}}/index/bwa/{{assembly}}".format(**config),
        algorithm=config["bwa_index_algo"]
    wrapper:
        "0.31.1/bio/bwa/index"


# TODO: maybe no wrapper, but make use of pipes/groups?
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#defining-groups-for-execution
# TODO merge params
rule bwa_mem:
    input:
        reads=expand("{fastq_dir}/trimmed/{{sample}}_trimmed.fastq.gz", **config),
        index=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
    output:
        expand("{result_dir}/mapped/{{sample}}-{{assembly}}.bam", **config)
    log:
        "logs/bwa_mem/{sample}-{assembly}.log"
    params:
        index=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}", **config),
        sort=config['bwa_mem_sort'],
        sort_order=config['bwa_mem_sort_order']
    priority: 1
    threads: 20
    wrapper:
        "file:../../wrappers/bwa/mem"


rule mark_duplicates:
    input:
        rules.bwa_mem.output
    output:
        bam=    expand("{result_dir}/dedup/{{sample}}-{{condition}}-{{project}}-{{assembly}}.bam", **config),
        metrics=expand("{result_dir}/dedup/{{sample}}-{{condition}}-{{project}}-{{assembly}}.metrics.txt", **config)
    log:
        "logs/mark_duplicates/{sample}-{condition}-{project}-{assembly}.log"
    params:
        config['duplicate_params']
    conda:
        "../envs/alignment.yaml"
    shell:
        "picard MarkDuplicates {params} INPUT={input} "
        "OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1"


rule samtools_stats:
    input:
        expand("{result_dir}/dedup/{{sample}}-{{condition}}-{{project}}-{{assembly}}.bam", **config)
    output:
        expand("{result_dir}/samtools_stats/{{sample}}-{{condition}}-{{project}}-{{assembly}}.txt", **config)
    log:
        "logs/samtools_stats/{sample}-{condition}-{project}-{assembly}.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"
