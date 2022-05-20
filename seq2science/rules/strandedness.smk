"""
all rules/logic related to the RNA-seq strandedness should be here.
"""
import os


rule infer_strandedness:
    """
    Use RSeqQC's infer_experiment.py to determine the strandedness of a sample
    """
    input:
        bam=FINAL_BAM,
        bai=FINAL_BAI,
        bed=rules.get_genome_annotation.output.bed,
    output:
        expand("{qc_dir}/strandedness/{{assembly}}-{{sample}}.strandedness.txt", **config),
    log:
        expand("{log_dir}/counts_matrix/{{assembly}}-{{sample}}.strandedness.log", **config),
    message: EXPLAIN["infer_strandedness"]
    params:
        config["min_mapping_quality"],
    conda:
        "../envs/gene_counts.yaml"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {input.bed} -q {params} 2> {log} | awk NF > {output} 
        """


def get_strandedness(report_file, fmt="htseq"):
    """
    Read RSeQC strandedness info into HTSeq/featurecount format.
    """
    if not os.path.exists(report_file):
        return "placeholder"

    cutoff = 0.6
    with open(report_file) as rf:
        fwd, rev = rf.read().splitlines()[2:4]

    fwd = float(fwd.split(": ")[1])
    if fwd > cutoff:
        return "yes" if fmt=="htseq" else 1
    rev = float(rev.split(": ")[1])
    if rev > cutoff:
        return "reverse" if fmt=="htseq" else 2
    return "no" if fmt=="htseq" else 0
