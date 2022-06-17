"""
all rules/logic related to motif scanning of peaks is found here.
"""
def get_motif2factors_input_genomes(wildcards):
    all_out = []
    if str(wildcards.assembly)[0:4] in ["GRCh", "GRCm", "hg19", "hg38", "mm10", "mm39"]:
        return all_out

    for assembly in config.get("motif2factors_database_references", []) + config.get("motif2factors_reference", []) + [wildcards.assembly]:
        all_out.append(expand(f"{{genome_dir}}/{assembly}/{assembly}.annotation.gtf", **config)[0])
        all_out.append(expand(f"{{genome_dir}}/{assembly}/{assembly}.fa", **config)[0])
    return all_out


rule motif2factors:
    """
    Create a gimme pfm/motif2factor file with the gimme motif2factors command.
    For human/mouse it just uses the default database, for other species it is based
    on TF orthologs.
    """
    input:
        genome=rules.get_genome.output,
        all_genomes=get_motif2factors_input_genomes
    output:
        expand("{result_dir}/gimme/{{assembly}}.{{gimme_database}}.pfm", **config),
    log:
        expand("{log_dir}/gimme/motif2factors/{{assembly}}-{{gimme_database}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/gimme/motif2factors/{{assembly}}-{{gimme_database}}.log", **config)[0],
    params:
        genomes_dir=config.get("genome_dir"),
        database=config.get("gimme_maelstrom_database"),
        motif2factors_reference=config.get("motif2factors_reference")
    threads: 24
    conda:
        "../envs/gimme.yaml"
    script:
        f"{config['rule_dir']}/../scripts/motif2factors.py"


rule gimme_maelstrom:
    """
    Gimme maelstrom is a method to infer differential motifs between two or more biological replicates.
    """
    input:
        genome=rules.get_genome.output,
        count_table=expand("{counts_dir}/{{peak_caller}}/{{assembly}}_log2_quantilenorm_biological_reps.tsv", **config),
        pfm=rules.motif2factors.output
    output:
        directory(expand("{result_dir}/gimme/maelstrom/{{assembly}}-{{gimme_database}}-{{peak_caller}}", **config)),
    params: config.get("gimme_maelstrom_params", "")
    log:
        expand("{log_dir}/gimme_maelstrom/{{assembly}}-{{gimme_database}}-{{peak_caller}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/gimme_maelstrom/{{assembly}}-{{gimme_database}}-{{peak_caller}}.log", **config)[0],
    resources:
        mem_gb=40,
    message: EXPLAIN["gimme_maelstrom"]
    conda:
        "../envs/gimme.yaml"
    threads: 24
    shell:
        ("cpulimit --include-children -l {threads}00 --\\" if config.get("cpulimit", True) else "") +
        """
        gimme maelstrom {input.count_table} {input.genome} {output} --pfmfile {input.pfm} --nthreads {threads} {params} > {log} 2>&1
        """
