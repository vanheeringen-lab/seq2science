rule infer_strandedness:
    """
    use RSeqQC's infer_experiment.py to determine strandedness af a sample
    """
    input:
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
        bai=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam.bai", **config),
        bed=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.bed", **config)
    output:
        expand("{qc_dir}/strandedness/{{assembly}}-{{sample}}.strandedness.txt", **config)
    log:
        expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.strandedness.log", **config),
    message: explain_rule("infer_strandedness")
    params:
        config["min_mapping_quality"]
    conda:
        "../envs/gene_counts.yaml"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} -q {params} 2> {log} | awk NF > {output} 
        """


def samples_to_infer(wildcards):
    """
    list all samples for which strandedness must be inferred
    """
    col = samples.technical_replicate if "technical_replicate" in samples else samples.index
    if config['ignore_strandedness'] or \
            ("strandedness" in samples and "nan" not in set(samples.strandedness)):
        files = []
    elif "strandedness" not in samples:
        files = [f"{{qc_dir}}/strandedness/{samples[col == sample].assembly[0]}{suffix}-{sample}.strandedness.txt" for sample in set(col)]
    else:
        files = []
        for sample in set(col):
            if samples[col == sample].strandedness not in ["yes", "forward", "reverse", "no"]:
                files.append(f"{{qc_dir}}/strandedness/{samples[col == sample].assembly[0]}{suffix}-{sample}.strandedness.txt")
    return expand(files, **config)


checkpoint strandedness_report:
    """
    combine samples.tsv & infer_strandedness results (call strandedness if >60% of reads explains a direction)
    """
    input:
        samples_to_infer
    output:
        expand("{qc_dir}/strandedness/inferred_strandedness.tsv", **config)
    params:
        reps=lambda wildcards, input: input  # help resolve changes in input files
    run:
        import pandas as pd

        def get_strand(sample):
            report_file = [f for f in input if f.endswith(f"-{sample}.strandedness.txt")][0]
            with open(report_file) as report:
                fail_val = fwd_val = 0
                for line in report:
                    if line.startswith("Fraction of reads failed"):
                        fail_val = float(line.strip().split(": ")[1])
                    elif line.startswith(("""Fraction of reads explained by "1++""",
                                          """Fraction of reads explained by "++""")):
                        fwd_val = float(line.strip().split(": ")[1])

            if fwd_val > 0.6:
                return "forward"
            elif 1 - (fwd_val + fail_val) > 0.6:
                return "reverse"
            else:
                return "no"

        strands = []
        method = []
        col = samples.technical_replicate if "technical_replicate" in samples else samples.index
        for sample in set(col):
            s = samples[col == sample].strandedness[0] if "strandedness" in samples else "nan"
            m = "user_specification"
            if config['ignore_strandedness']:
                s = "no"
                m = "ignored"
            elif s == "nan":
                s = get_strand(sample)
                m = "inferred"
            strands.append(s)
            method.append(m)

        strandedness = pd.DataFrame({"sample": list(set(col)), "strandedness": strands, "determined_by": method}, dtype='str')
        strandedness.set_index('sample', inplace=True)
        strandedness.to_csv(output[0], sep="\t")


def _strandedness_report(wildcards):
    strandreport = checkpoints.strandedness_report.get().output
    return strandreport[0]


def strandedness_to_quant(wildcards, tool):
    """
    translate strandedness to quantifiers nomenclature
    """
    out = {
        "htseq": ["no", "yes", "reverse"],
        "featurecounts": ["0", "1", "2"],
        "dexseq": ["no", "yes", "reverse"],
    }

    strandedness = pd.read_csv(_strandedness_report(wildcards), sep='\t', dtype='str', index_col=0)
    if wildcards.sample not in strandedness.index:
        return "new samples added, start rerun"
    s = strandedness[strandedness.index == wildcards.sample].strandedness[0]
    n = 1 if s in ["yes", "forward"] else (2 if s == "reverse" else 0)
    return out[tool][n]
