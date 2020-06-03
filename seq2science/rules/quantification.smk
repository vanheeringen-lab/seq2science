if config['quantifier'] == 'star':
    # star_align also generates bams, required for the trackhub
    if config['create_trackhub']:
        ruleorder: star_align > star_quant
    else:
        ruleorder: star_quant > star_align


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
            directory(expand("{log_dir}/{quantifier}_quant/{{assembly}}-{{sample}}", **config))
        benchmark:
            expand("{benchmark_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f' {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f' {input.reads[0]} {input.reads[1]}',
            params=config['quantify']
        threads: 8
        resources:
            mem_gb=30
        conda:
            "../envs/star.yaml"
        shell:
            """
            trap "find {log} -type f ! -name Log* -exec rm {{}} \;" EXIT
            mkdir -p {log}
            mkdir -p {output.dir}                
            
            STAR --genomeDir {input.index} --readFilesIn {params.input} --readFilesCommand gunzip -c \
            --quantMode GeneCounts --outSAMtype None \
            --outFileNamePrefix {log}/ --outTmpDir {output.dir}/STARtmp \
            --runThreadN {threads} {params.params} > {log}/Log.std_stderr.out 2>&1
            
            # move all non-log files to output directory (this way the log files are kept on error)
            find {log} -type f ! -name Log* -exec mv {{}} {output.dir} \;
            """


elif config['quantifier'] == 'salmon':
    rule salmon_decoy_aware_index:
        """
        Make a decoy aware transcript index for Salmon.
        """
        input:
            transcripts=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
            decoy_transcripts=expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config),
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{quantifier}_decoy_aware", **config))
        log:
            expand("{log_dir}/{quantifier}_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{quantifier}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config['quantifier_index']
        threads: 40
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon index -t {input.transcripts} --decoys {input.decoy_transcripts} -i {output} {params} \
            --threads {threads} > {log} 2>&1
            """

    rule salmon_index:
        """
        Make a transcriptomic index for Salmon.
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
            config['quantifier_index']
        threads: 40
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon index -t {input} -i {output} {params} --threads {threads} > {log} 2>&1
            """

    rule salmon_quant:
        """
        Align reads against a transcriptome (index) with Salmon (mapping-based mode) and output a quantification file per sample.
        """
        input:
            reads=get_reads,
            index=get_index
        output:
            dir=directory(expand("{result_dir}/{quantifier}/{{assembly}}-{{sample}}", **config)),
        log:
            expand("{log_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-r {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['quantify']
        threads: 20
        resources:
            mem_gb=8
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon quant -i {input.index} -l A {params.input} {params.params} -o {output.dir} \
            --threads {threads} > {log} 2>&1
            """
