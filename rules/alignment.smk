def get_reads(wildcards):
    if config.get('combine_replicates', '') == 'merge':
        if config['layout'].get(wildcards.sample, False) == "SINGLE":
            return expand("{result_dir}/{trimmed_dir}/merged/{{sample}}_trimmed.{fqsuffix}.gz", **config)
        return sorted(expand("{result_dir}/{trimmed_dir}/merged/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))
    else:
        if config['layout'].get(wildcards.sample, False) == "SINGLE":
            return expand("{result_dir}/{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config)
        return sorted(expand("{result_dir}/{trimmed_dir}/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))

def get_alignment_pipes():
    pipes = set()
    if config.get('peak_caller', False):
        if 'macs2' in config['peak_caller'] or 'hmmratac' in config['peak_caller']:
            pipes.add(pipe(expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.samtools.pipe", **config)[0]))
        if 'genrich' in config['peak_caller']:
            pipes.add(pipe(expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.sambamba.pipe", **config)[0]))
    else:
        pipes.add(pipe(expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.{bam_sorter}.pipe", **config)[0]))

    return pipes


if config['aligner'] == 'bowtie2':
    rule bowtie2_index:
        """
        Make a genome index for bowtie2.
        """
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
        params: config['index']
        shell:
            "bowtie2-build {params} --threads {threads} {input} {output}/{wildcards.assembly} > {log} 2>&1"


    rule bowtie2_align:
        """
        Align reads against a genome (index) with bowtie2, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bowtie2/", **config)
        output:
            get_alignment_pipes()
        log:
            expand("{log_dir}/bowtie2_align/{{sample}}-{{assembly}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bowtie2_align/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['params']
        threads: 20
        conda:
            "../envs/bowtie2.yaml"
        shell:
            """
            bowtie2 {params.params} --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'bwa':
    config['bwaindex_types'] = ['amb', 'ann', 'bwt', 'pac', 'sa']

    rule bwa_index:
        """
        Make a genome index for bwa.
        """
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
            params=config['index']
        conda:
            "../envs/bwa.yaml"
        shell:
            "bwa index -p {params.prefix} {params.params} {input} > {log} 2>&1"


    rule bwa_mem:
        """
        Align reads against a genome (index) with bwa, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
        output:
            get_alignment_pipes()
        log:
            expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bwa_mem/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            index_dir=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}", **config),
            params=config['align']
        threads: 20
        conda:
            "../envs/bwa.yaml"
        shell:
            """
            bwa mem {params.params} -t {threads} {params.index_dir} {input.reads} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'hisat2':
    rule hisat2_index:
        """
        Make a genome index for hisat2.
        """
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
        params: config['index']
        shell:
            "hisat2-build {params} -p {threads} {input} {output}/{wildcards.assembly} > {log} 2>&1"


    rule hisat2_align:
        """
        Align reads against a genome (index) with hisat2, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/hisat2/", **config)
        output:
            get_alignment_pipes()
        log:
            expand("{log_dir}/hisat2_align/{{sample}}-{{assembly}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/hisat2_align/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['align']
        threads: 20
        conda:
            "../envs/hisat2.yaml"
        shell:
            """
            hisat2 {params.params} --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'salmon':
    rule salmon_index:
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{aligner}", **config))
        log:
            expand("{log_dir}/{aligner}_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{aligner}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config['index']
        threads: 4
        conda:
            "../envs/salmon.yaml"
        shell:
            "salmon index -t {input} -i {output} {params} --threads {threads} &> {log}"


    rule salmon_quant:
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/{aligner}", **config)
        output:
            dir=directory(expand("{result_dir}/{aligner}/{{assembly}}/{{sample}}", **config)), #this could become a temp() directory, but quant.sf files are useful for other (currently unsupported) analyses
            pipe=get_alignment_pipes()
        log:
            expand("{log_dir}/{aligner}_align/{{sample}}-{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{aligner}_align/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-r {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['align']
        threads: 20
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon quant -i {input.index} -l A {params.input} {params.params} -o {output.dir} \
            --threads $(expr 4 * {threads} / 5) --writeMappings 2> {log} | \
            samtools view -b - -@ $(expr {threads} / 5) | tee {output.pipe} 1> /dev/null 2>> {log}
            """


rule sambamba_sort:
    """
    Sort the result of alignment with the sambamba sorter.
    """
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
    """
    Sort the result of alignment with the samtools sorter.
    """
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
    """
    Create an index of a bam file which can be used for e.g. visualization.
    """
    input:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{{bam_sorter}}-{{sorting}}.bam", **config)
    output:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{{bam_sorter}}-{{sorting}}.bai", **config)
    log:
        expand("{log_dir}/samtools_index/{{sample}}-{{assembly}}-{{bam_sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/samtools_index/{{sample}}-{{assembly}}-{{bam_sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        config['samtools_index']
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {params} {input} {output}
        """


rule mark_duplicates:
    """
    Mark (but keep) all duplicate reads in a bam file with picard MarkDuplicates
    """
    input:
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.{{sorter}}-{{sorting}}.bam", **config)
    output:
        bam=    expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{{sorter}}-{{sorting}}.bam", **config),
        metrics=expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.{{sorter}}-{{sorting}}.metrics.txt", **config)
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
