def get_reads(wildcards):
    if config['layout'][wildcards.sample] == "SINGLE":
        return expand("{result_dir}/{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config)
    return sorted(expand("{result_dir}/{trimmed_dir}/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))


if config['aligner'] == 'bowtie2':
    rule bowtie2_index:
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/bowtie2/", **config))
        log:
            expand("{log_dir}/bowtie2_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/bowtie2_index/{{assembly}}.benchmark.txt", **config)[0]
        threads: 4
        conda:
            "../envs/bowtie2.yaml"
        shell:
            "bowtie2-build --threads {threads} {input} {output}/{wildcards.assembly} > {log} 2>&1"


    rule bowtie2_align:
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bowtie2/", **config)
        output:
            pipe(expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bampipe", **config))
        log:
            expand("{log_dir}/bowtie2_align/{{sample}}-{{assembly}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bowtie2_align/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}'
        threads: 20
        conda:
            "../envs/bowtie2.yaml"
        shell:
            """
            bowtie2 --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 1> {output} 2> {log}
            """

elif config['aligner'] == 'bwa':
    config['bwaindex_types'] = ['amb', 'ann', 'bwt', 'pac', 'sa']

    rule bwa_index:
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
        log:
            expand("{log_dir}/bwa_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/bwa_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            prefix="{genome_dir}/{{assembly}}/index/bwa/{{assembly}}".format(**config),
            algorithm=config["bwa_algo"]
        conda:
            "../envs/bwa.yaml"
        shell:
            "bwa index -p {params.prefix} -a {params.algorithm} {input} > {log} 2>&1"

    def get_pipes():
        pipes = set()
        if config.get('peak_caller', False):
            if 'macs2' in config['peak_caller']:
                pipes.add(pipe(expand("{result_dir}/{bwa_dir}/{{sample}}-{{assembly}}.samtools.pipe", **config)[0]))
            if 'genrich' in config['peak_caller']:
                pipes.add(pipe(expand("{result_dir}/{bwa_dir}/{{sample}}-{{assembly}}.sambamba.pipe", **config)[0]))
            if 'hmmratac' in config['peak_caller']:
                pipes.add(pipe(expand("{result_dir}/{bwa_dir}/{{sample}}-{{assembly}}.samtools.pipe", **config)[0]))
        else:
            pipes.add(pipe(expand("{result_dir}/{bwa_dir}/{{sample}}-{{assembly}}.{bam_sorter}.pipe", **config)[0]))

        return pipes

    rule bwa_mem:
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
        output:
            get_pipes()
        log:
            expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bwa_mem/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            index_dir=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}", **config)
        threads: 20
        conda:
            "../envs/bwa.yaml"
        shell:
            """
            bwa mem -t {threads} {params.index_dir} {input.reads} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """

elif config['aligner'] == 'hisat2':
    rule hisat2_index:
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/hisat2/", **config))
        log:
            expand("{log_dir}/hisat2_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/hisat2_index/{{assembly}}.benchmark.txt", **config)[0]
        threads: 4
        conda:
            "../envs/hisat2.yaml"
        shell:
            "hisat2-build -p {threads} {input} {output}/{wildcards.assembly} > {log} 2>&1"


    rule hisat2_align:
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/hisat2/", **config)
        output:
            pipe(expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bampipe", **config))
        log:
            expand("{log_dir}/hisat2_align/{{sample}}-{{assembly}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/hisat2_align/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}'
        threads: 20
        conda:
            "../envs/hisat2.yaml"
        shell:
            """
            hisat2 --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 1> {output} 2> {log}
            """


rule sambamba_sort:
    input:
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.sambamba.pipe", **config)
    output:
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.sambamba-{{sorting}}.bam", **config)
    log:
        expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}-sambamba_{{sorting}}.log", **config)
    group: 'alignment'
    benchmark:
        expand("{benchmark_dir}/sambamba_sort/{{sample}}-{{assembly}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        lambda wildcards: "-n" if wildcards.sorting == 'queryname' else '',
    threads: 4
    conda:
        "../envs/sambamba.yaml"
    shell:
        """
        sambamba view --nthreads {threads} -S -f bam  {input[0]} -o /dev/stdout  2> {log} |
        sambamba sort --nthreads {threads} {params}   /dev/stdin -o {output[0]}  2> {log}
        """


rule samtools_sort:
    input:
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.samtools.pipe", **config)
    output:
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.samtools-{{sorting}}.bam", **config)
    log:
        expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}-samtools_{{sorting}}.log", **config)
    group: 'alignment'
    benchmark:
        expand("{benchmark_dir}/samtools_sort/{{sample}}-{{assembly}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        order=lambda wildcards: "-n" if wildcards.sorting == 'queryname' else '',
        threads=lambda wildcards, input, output, threads: threads - 1
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        trap \"rm -f {output}*\" INT;
        samtools sort -@ {params.threads} {params.order} {input} -o {output}  2> {log}
        """


rule samtools_index:
    input:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{bam_sorter}-{{sorting}}.bam", **config)
    output:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{bam_sorter}-{{sorting}}.bai", **config)
    log:
        expand("{log_dir}/samtools_index/{{sample}}-{{assembly}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/samtools_index/{{sample}}-{{assembly}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        config['samtools_index']
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {params} {input} {output}
        """


rule mark_duplicates:
    input:
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.{{sorter}}-{{sorting}}.bam", **config)
    output:
        bam=    expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{{sorter}}-{{sorting}}.bam", **config),
        metrics=expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}-{{sorter}}-{{sorting}}.metrics.txt", **config)
    log:
        expand("{log_dir}/mark_duplicates/{{sample}}-{{assembly}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/mark_duplicates/{{sample}}-{{assembly}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        config['markduplicates']
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates {params} INPUT={input} "
        "OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1"
