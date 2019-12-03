# STAR can produce bams and gene counts at the same time. This will be done by the rules in alignment,smk
if config['quantifier'] == 'star' and not config['run_alignment']:
    # rule star_index can be found in alignment.smk
    rule star_quant:
        """
        Quantify reads against a genome and transcriptome (index) with STAR and output a counts table per sample.
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/{quantifier}", **config)
        output:
            dir=directory(expand("{result_dir}/{quantifier}/{{assembly}}-{{sample}}", **config)),
        log:
            expand("{log_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f' {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f' {input.reads[0]} {input.reads[1]}',
            params=config['quantifier_flags']
        threads: 1
        resources:
            mem_gb=30
        conda:
            "../envs/star.yaml"
        shell:
            """
            mkdir -p {output.dir}
            
            STAR --genomeDir {input.index} --readFilesIn {params.input} --quantMode GeneCounts \
            --outFileNamePrefix {output.dir}/ --runThreadN {threads} {params.params} \
            --outSAMtype None &> {log}
            """


elif config['quantifier'] == 'salmon':
    rule salmon_index:
        """
        Make a transcript index for Salmon.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{quantifier}", **config))
        log:
            expand("{log_dir}/{quantifier}_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{quantifier}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config['quantifier_index_flags']
        threads: 40
        conda:
            "../envs/salmon.yaml"
        shell:
            "salmon index -t {input} -i {output} {params} --threads {threads} &> {log}"


    rule salmon_quant:
        """
        Align reads against a transcriptome (index) with Salmon (mapping-based mode) and output a quantification file per sample.
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/{quantifier}", **config)
        output:
            dir=directory(expand("{result_dir}/{quantifier}/{{assembly}}-{{sample}}", **config)),
        log:
            expand("{log_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-r {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['quantifier_flags']
        threads: 20
        resources:
            mem_gb=8
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon quant -i {input.index} -l A {params.input} {params.params} -o {output.dir} \
            --threads $(( 4 * {threads} / 5)) 2> {log}
            """
