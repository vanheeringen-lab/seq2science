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


    rule bwa_mem:
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
        output:
            pipe(expand("{result_dir}/{bwa_dir}/{{sample}}-{{assembly}}.bampipe", **config))
        log:
            expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/bwa_mem/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            index_dir=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}", **config)
        threads: 20
        conda:
            "../envs/bwa.yaml"
        shell:
            """
            bwa mem -t {threads} {params.index_dir} {input.reads} 1> {output} 2> {log}
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


elif config['aligner'] == 'salmon':
    rule salmon_index:
        input:
            #expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
            expand("{genome_dir}/{{assembly}}/transcriptome.fa", **config) #todo: expand genomepy to download transcriptomes
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{aligner}/", **config))
        log:
            expand("{log_dir}/{aligner}_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{aligner}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config['salmon_ind']
        threads: 4
        conda:
            "../envs/salmon.yaml"
        shell:
            "salmon index -t {input} -i {output} {params} --threads {threads} > {log} 2&1"

    rule salmon_quant:
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/{aligner}/", **config)
        output:
            pipe(expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bampipe", **config))
        log:
            expand("{log_dir}/{aligner}_align/{{sample}}-{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{aligner}_align/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-r {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            flags=config['salmon_aln']
        threads: 21
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon quant -i {input.index} -l A {params.input} {params.flags} --threads {threads}/3 2> {log} --writeMappings | samtools view --threads {threads}/3 -Sb - | samtools sort --threads {threads}/3 -T sort.tmp -o - > {output}
            """

            # """
            # salmon quant -i {input.index} {params.input} 2> {log} \
            # # output stdout in SAM format, and convert this to BAM.
            # # requires mapping-based mode with a quasi-index
            # --writeMappings | samtools view -Sb - | samtools sort -T sort.tmp -o - > {output}
            # """
            #
            # "salmon quant -i {params.salmonindex} -l A -1 {input.fastq1} -2 {input.fastq2} "
		    # "-o {params.salmondir}/{wildcards.sample} --seqBias --gcBias "
		    # "--fldMean {params.fldMean} --fldSD {params.fldSD} -p {threads}"
		    # #"--posBias"


if config['bam_sorter'] == 'sambamba':
    rule sambamba_sort:
        input:
            expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bampipe", **config)
        output:
            expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bam", **config)
        log:
            expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}-sambamba_sort.log", **config)
        benchmark:
            expand("{benchmark_dir}/sambamba_sort/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            "-n" if config['bam_sort_order'] == 'queryname' else ''
        threads: 4
        conda:
            "../envs/sambamba.yaml"
        shell:
            """
            sambamba view --nthreads {threads} -S -f bam  {input[0]} -o /dev/stdout  2> {log} |
            sambamba sort --nthreads {threads} {params}   /dev/stdin -o {output[0]}  2> {log}
            """


elif config['bam_sorter'] == 'samtools':
    rule samtools_sort:
        input:
            expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bampipe", **config)
        output:
            expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bam", **config)
        log:
            expand("{log_dir}/bwa_mem/{{sample}}-{{assembly}}-samtools_sort.log", **config)
        benchmark:
            expand("{benchmark_dir}/samtools_sort/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
        params:
            order="-n" if config['bam_sort_order'] == 'queryname' else '',
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
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bam", **config)
    output:
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bai", **config)
    log:
        expand("{log_dir}/samtools_index/{{sample}}-{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/samtools_index/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
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
        expand("{result_dir}/{aligner}/{{sample}}-{{assembly}}.bam", **config)
    output:
        bam=    expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.bam", **config),
        metrics=expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.metrics.txt", **config)
    log:
        expand("{log_dir}/mark_duplicates/{{sample}}-{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/mark_duplicates/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
    params:
        config['markduplicates']
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates {params} INPUT={input} "
        "OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1"
