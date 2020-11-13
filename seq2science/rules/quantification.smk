if config["quantifier"] == "salmon":

    rule get_transcripts:
        """
        Generate transcripts.fasta using gffread.
    
        Requires genome.fa and annotation.gtf (with matching chromosome/scaffold names)
        """
        input:
            fa=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        output:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
        log:
            expand("{log_dir}/get_genome/{{assembly}}.transcripts.log", **config),
        benchmark:
            expand("{benchmark_dir}/get_genome/{{assembly}}.transcripts.benchmark.txt", **config)[0]
        conda:
            "../envs/salmon.yaml"
        priority: 1
        shell:
            "gffread -w {output} -g {input.fa} {input.gtf} >> {log} 2>&1"


    rule decoy_transcripts:
        """
        Generate decoy_transcripts.txt for Salmon indexing  
    
        script source: https://github.com/COMBINE-lab/SalmonTools
        """
        input:
            genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            transcripts=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
        output:
            gentrome=expand("{genome_dir}/{{assembly}}/decoy_transcripts/gentrome.fa", **config),
            decoys=expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config),
        params:
            script=f"{config['rule_dir']}/../scripts/generateDecoyTranscriptome.sh",
        log:
            expand("{log_dir}/get_genome/{{assembly}}.decoy_transcripts.log", **config),
        message: explain_rule("decoy_transcripts")
        benchmark:
            expand("{benchmark_dir}/get_genome/{{assembly}}.decoy_transcripts.benchmark.txt", **config)[0]
        threads: 40
        resources:
            mem_gb=65,
        conda:
            "../envs/decoy.yaml"
        priority: 1
        shell:
            ("cpulimit --include-children -l {threads}00 -- " if config. get("cpulimit", True) else" ")+
            "sh {params.script} -j {threads} -g {input.genome} -a {input.gtf} -t {input.transcripts} -o $(dirname {output[0]}) > {log} 2>&1"

    rule salmon_decoy_aware_index:
        """
        Generate a decoy aware transcriptome index for Salmon.
        """
        input:
            gentrome=expand("{genome_dir}/{{assembly}}/decoy_transcripts/gentrome.fa", **config),
            decoys=expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config),
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{quantifier}_decoy_aware", **config)),
        log:
            expand("{log_dir}/{quantifier}_index/{{assembly}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/{quantifier}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config["quantifier_index"],
        priority: 1
        threads: 10
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon index -t {input.gentrome} -d {input.decoys} -i {output} {params} \
            --threads {threads} > {log} 2>&1
            """

    rule salmon_index:
        """
        Generate a transcriptome index for Salmon.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{quantifier}", **config)),
        log:
            expand("{log_dir}/{quantifier}_index/{{assembly}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/{quantifier}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config["quantifier_index"],
        priority: 1
        threads: 10
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
            index=get_salmon_index,
        output:
            dir=directory(expand("{result_dir}/{quantifier}/{{assembly}}-{{sample}}", **config)),
        log:
            expand("{log_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.log", **config),
        message: explain_rule("salmon_quant")
        benchmark:
            expand("{benchmark_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=(
                lambda wildcards, input: ["-r", input.reads]
                if sampledict[wildcards.sample]["layout"] == "SINGLE"
                else ["-1", input.reads[0], "-2", input.reads[1]]
            ),
            params=config["quantifier_flags"],
        threads: 12
        resources:
            mem_gb=8,
        conda:
            "../envs/salmon.yaml"
        shell:
            """
            salmon quant -i {input.index} -l A {params.input} {params.params} -o {output.dir} \
            --threads {threads} > {log} 2>&1
            """


elif config["quantifier"] == "htseq":

    rule htseq_count:
        """
        summarize reads to gene level. Outputs a counts table per bam file.
        """
        input:
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            required=_strandedness_report,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv", **config),
        params:
            strandedness=lambda wildcards: strandedness_to_quant(wildcards, "htseq"),
            user_flags=config["htseq_flags"]
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.counts.log", **config),
        message: explain_rule("htseq_count")
        threads: 1
        conda:
            "../envs/gene_counts.yaml"
        shell:
             """
             htseq-count {input.bam} {input.gtf} -r pos -s {params.strandedness} {params.user_flags} -n {threads} -c {output} > {log} 2>&1
             """


elif config["quantifier"] == "featurecounts":

    rule featurecounts:
        """
        summarize reads to gene level. Outputs a counts table per bam file.
        """
        input:
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            required=_strandedness_report,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv", **config),
        params:
            strandedness=lambda wildcards: strandedness_to_quant(wildcards, "featurecounts"),
            endedness=lambda wildcards: "" if sampledict[wildcards.sample]["layout"] == 'SINGLE' else "-p",
            user_flags=config["featurecounts_flags"],
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.counts.log", **config),
        message: explain_rule("featurecounts_rna")
        threads: 1
        conda:
            "../envs/gene_counts.yaml"
        shell:
            """
            featureCounts -a {input.gtf} {input.bam} {params.endedness} -s {params.strandedness} {params.user_flags} -T {threads} -o {output} > {log} 2>&1
            """


if config.get("dexseq"):

    rule prepare_DEXseq_annotation:
        """
        generate a DEXseq annotation.gff from the annotation.gtf
        """
        input:
             expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        output:
             expand("{genome_dir}/{{assembly}}/{{assembly}}.DEXseq_annotation.gff", **config),
        log:
             expand("{log_dir}/counts_matrix/{{assembly}}.prepare_DEXseq_annotation.log", **config),
        conda:
             "../envs/dexseq.yaml"
        shell:
             """
             current_conda_env=$(conda env list | grep \* | cut -d "*" -f2-)
             DEXseq_path=${{current_conda_env}}/lib/R/library/DEXSeq/python_scripts
             
             python ${{DEXseq_path}}/dexseq_prepare_annotation.py {input} {output} > {log} 2>&1
             """

    rule dexseq_count:
        """
        count exon usage
        """
        input:
            bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
            bai=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam.bai", **config),
            gff=expand("{genome_dir}/{{assembly}}/{{assembly}}.DEXseq_annotation.gff", **config),
            required=_strandedness_report,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.DEXSeq_counts.tsv", **config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.DEXseq_counts.log", **config),
        message: explain_rule("dexseq")
        params:
            strandedness=lambda wildcards: strandedness_to_quant(wildcards, "dexseq"),
            endedness=lambda wildcards: "" if sampledict[wildcards.sample]["layout"] == 'SINGLE' else "-p yes",
        conda:
             "../envs/dexseq.yaml"
        shell:
             """
             current_conda_env=$(conda env list | grep \* | cut -d "*" -f2-)
             DEXseq_path=${{current_conda_env}}/lib/R/library/DEXSeq/python_scripts
             
             python ${{DEXseq_path}}/dexseq_count.py -f bam -r pos {params.endedness} -s {params.strandedness} {input.gff} {input.bam} {output} > {log} 2>&1
             """
