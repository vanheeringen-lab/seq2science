"""
all rules/logic related to the RNA-seq strandedness should be here.
"""


rule infer_strandedness:
    """
    use RSeqQC's infer_experiment.py to determine strandedness af a sample
    """
    input:
        bam=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config),
        bai=expand("{final_bam_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam.bai", **config),
        bed=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.bed", **config),
    output:
        expand("{qc_dir}/strandedness/{{assembly}}-{{sample}}.strandedness.txt", **config),
    log:
        expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.strandedness.log", **config),
    message:
        explain_rule("infer_strandedness")
    params:
        config["min_mapping_quality"],
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
    col = samples.technical_replicates if "technical_replicates" in samples else samples.index
    if config["ignore_strandedness"] or ("strandedness" in samples and "nan" not in set(samples.strandedness)):
        files = []
    elif "strandedness" not in samples:
        files = [
            f"{{qc_dir}}/strandedness/{samples[col == sample].assembly[0]}{suffix}-{sample}.strandedness.txt"
            for sample in set(col)
        ]
    else:
        files = []
        for sample in set(col):
            if samples[col == sample].strandedness not in ["yes", "forward", "reverse", "no"]:
                files.append(
                    f"{{qc_dir}}/strandedness/{samples[col == sample].assembly[0]}{suffix}-{sample}.strandedness.txt"
                )
    return list(sorted(expand(files, **config)))


rule strandedness_report:
    """
    combine samples.tsv & infer_strandedness results (call strandedness if >60% of reads explains a direction)
    """
    input:
        samples_to_infer,
    output:
        expand("{qc_dir}/strandedness/inferred_strandedness.tsv", **config),
    params:
        input=lambda wildcards, input: input,  # help resolve changes in input files
        reps=treps,
    run:
        import pandas as pd


        def get_strand(sample):
            report_file = [f for f in input if f.endswith(f"-{sample}.strandedness.txt")][0]
            with open(report_file) as report:
                fail_val = fwd_val = 0
                for line in report:
                    if line.startswith("Fraction of reads failed"):
                        fail_val = float(line.strip().split(": ")[1])
                    elif line.startswith(
                        ("""Fraction of reads explained by "1++""", """Fraction of reads explained by "++""")
                    ):
                        fwd_val = float(line.strip().split(": ")[1])

            if fwd_val > 0.6:
                return "yes"
            elif 1 - (fwd_val + fail_val) > 0.6:
                return "reverse"
            else:
                return "no"


        strands = []
        method = []
        col = samples.technical_replicates if "technical_replicates" in samples else samples.index
        for sample in set(col):
            s = samples[col == sample].strandedness[0] if "strandedness" in samples else "nan"
            m = "user_specification"
            if config["ignore_strandedness"]:
                s = "no"
                m = "ignored"
            elif s == "nan":
                s = get_strand(sample)
                m = "inferred"
            strands.append(s)
            method.append(m)

        strandedness = pd.DataFrame(
            {"sample": list(set(col)), "strandedness": strands, "determined_by": method}, dtype="str"
        )
        strandedness.set_index("sample", inplace=True)
        strandedness.to_csv(output[0], sep="\t")
