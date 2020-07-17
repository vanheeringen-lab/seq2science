def get_strandedness(wildcards):
    out = {
        "htseq": ["no", "yes", "reverse"],
        "featurecounts": ["0", "1", "2"]
    }

    if "strandedness" not in samples:
        n = 0
    else:
        col = "replicate" if "replicate" in samples else "sample"
        s = samples[samples[col] == wildcards.sample].strandedness[0]

        if s in ["yes", "forward"]:
            n=1
        elif s == "reverse":
            n=2
        else:
            n=0
    return out[config["quantifier"]][n]


if config["quantifier"] == "salmon":


    rule salmon_decoy_aware_index:
        """
        Generate a decoy aware transcriptome index for Salmon.
        """
        input:
            transcripts=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
            decoy_transcripts=expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config),
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{quantifier}_decoy_aware", **config)),
        log:
            expand("{log_dir}/{quantifier}_index/{{assembly}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/{quantifier}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config["quantifier_index"],
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
            index=get_index,
        output:
            dir=directory(expand("{result_dir}/{quantifier}/{{assembly}}-{{sample}}", **config)),
        log:
            expand("{log_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=(
                lambda wildcards, input: ["-r", input.reads]
                if config["layout"][wildcards.sample] == "SINGLE"
                else ["-1", input.reads[0], "-2", input.reads[1]]
            ),
            params=config["quantify"],
        threads: 20
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
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv", **config),
        params:
            strandedness=get_strandedness,
            user_flags=config["htseq_flags"]
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.counts.log", **config),
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
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv", **config),
        params:
            strandedness=get_strandedness,
            endedness=lambda wildcards: "" if config['layout'][wildcards.sample] == 'SINGLE' else "-p",
            user_flags=config["featurecounts_flags"],
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.counts.log", **config),
        threads: 1
        conda:
            "../envs/gene_counts.yaml"
        shell:
            """
            featureCounts -a {input.gtf} {input.bam} {params.endedness} -s {params.strandedness} {params.user_flags} -T {threads} -o {output} > {log} 2>&1
            """
