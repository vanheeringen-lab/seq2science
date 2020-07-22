def get_strandedness(wildcards):
    out = {
        "htseq": ["no", "yes", "reverse"],
        "featurecounts": ["0", "1", "2"]
    }
    strandedness = pd.read_csv(get_strand_table(wildcards), sep='\t', dtype='str', index_col=0)

    s = strandedness[strandedness.index == wildcards.sample].strandedness[0]
    if s in ["yes", "forward"]:
        n=1
    elif s == "reverse":
        n=2
    else:
        n=0
    # if "strandedness" not in samples:
    #     n = 0
    # else:
    #     col = samples.replicate if "replicate" in samples else samples.index
    #     s = samples[col == wildcards.sample].strandedness[0]
    #
    #     if s in ["yes", "forward"]:
    #         n=1
    #     elif s == "reverse":
    #         n=2
    #     else:
    #         n=0
    return out[config["quantifier"]][n]

rule infer_sample_strandedness:
    """
    use RSeqQC infer_experiment.py to determine strandedness if none was provided in the sample.tsv
    """
    input:
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
        bed=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.bed", **config)
    output:
        temp(expand("{counts_dir}/{{assembly}}-{{sample}}.strandedness.txt", **config))
    log:
        expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.strandedness.log", **config),
    params:
        config["min_mapping_quality"]
    conda:
        "../envs/gene_counts.yaml"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} -q {params} 1> {output} 2> {log}
        """

# infer_experiment.py \
# -i /bank/experiments/2020-07/colin/final_bam/hg19-external-CL-RNA10-P4-plusTA-18410.samtools-coordinate.bam \
# -r /home/siebrenf/.local/share/genomes/hg19/hg19.annotation.bed \
# -q 44

def strandedness_reports(wildcards):
    """
    list all samples for which strandedness must be inferred
    """
    col = samples.replicate if "replicate" in samples else samples.index

    if "strandedness" not in samples:
        files = [f"{{counts_dir}}/{samples[col == sample].assembly[0]}-{sample}.strandedness.txt" for sample in col]
    else:
        files = []
        for sample in samples:
            if samples[col == sample].strandedness not in ["yes", "forward", "reverse", "no"]:
                files.append(f"{{counts_dir}}/{samples[col == sample].assembly[0]}-{sample}.strandedness.txt")
    return expand(files, **config)


checkpoint infer_strandedness:
    input:
        strandedness_reports
    output:
        expand("{counts_dir}/inferred_strandedness.tsv", **config)
    run:
        import pandas as pd

        # samples = pd.read_csv("/bank/experiments/2020-07/colin/samples.tsv", sep='\t', dtype='str', comment='#')
        # samples = samples.set_index('sample')

        # populate the table with the samples.tsv
        col = samples.replicate if "replicate" in samples else samples.index
        strandedness = pd.DataFrame({"sample": list(col), "strandedness": samples.strandedness if "strandedness" in samples else "nan", "obtained": "from_samples.tsv"}, dtype=str)
        strandedness.set_index('sample', inplace=True)

        if "nan" in list(strandedness.strandedness):
            def get_strand(sample):
                report_file = [f for f in input if f.endswith(f"-{sample}.strandedness.txt")][0]
                with open(report_file) as report:
                    fail_val = fwd_val = 0
                    for line in report:
                        if line.startswith("Fraction of reads failed"):
                            fail_val = float(line.strip().split(": ")[1])
                        elif line.startswith("""Fraction of reads explained by "1++"""):
                            fwd_val = float(line.strip().split(": ")[1])
                        elif line.startswith("""Fraction of reads explained by "++"""):
                            fwd_val = float(line.strip().split(": ")[1])

                if fwd_val > 0.6:
                    return "forward"
                elif 1 - (fwd_val + fail_val) > 0.6:
                    return "reverse"
                else:
                    return "no"

            strands = []
            inf = []
            for sample in strandedness.index:
                s = strandedness[strandedness.index == sample].strandedness[0]
                i = "from_samples.tsv"
                if s == "nan":
                    s = get_strand(sample)
                    i = "inferred"
                strands.append(s)
                inf.append(i)

            strandedness = pd.DataFrame({"sample": list(col), "strandedness": strands, "obtained": inf}, dtype=str)
            strandedness.set_index('sample', inplace=True)

        strandedness.to_csv(output[0], sep="\t")

def get_strand_table(wildcards):
    return checkpoints.infer_strandedness.rule.output[0]

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
            expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config),
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
            "sh {params.script} -j {threads} -g {input.genome} -a {input.gtf} -t {input.transcripts} -o $(dirname {output}) > {log} 2>&1"

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
        message: explain_rule("salmon_quant")
        benchmark:
            expand("{benchmark_dir}/{quantifier}_quant/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=(
                lambda wildcards, input: ["-r", input.reads]
                if config["layout"][wildcards.sample] == "SINGLE"
                else ["-1", input.reads[0], "-2", input.reads[1]]
            ),
            params=config["quantifier_flags"],
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
            test=get_strand_table,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv", **config),
        params:
            strandedness=get_strandedness,
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
            test=get_strand_table,
        output:
            expand("{counts_dir}/{{assembly}}-{{sample}}.counts.tsv", **config),
        params:
            strandedness=get_strandedness,
            endedness=lambda wildcards: "" if config['layout'][wildcards.sample] == 'SINGLE' else "-p",
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
