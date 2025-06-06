"""
all rules/logic related to downloading public assemblies and preprocessing of assemblies should be here.
"""


localrules:
    extend_genome,
    get_genome_support_files,


rule get_genome:
    """
    Download a genome through genomepy.
    """
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
    log:
        expand("{log_dir}/get_genome/{{assembly}}.genome.log", **config),
    message: EXPLAIN["get_genome"]
    params:
        providers=PROVIDERS,
        provider=config.get("provider"),
        genomes_dir=config["genome_dir"],
    wildcard_constraints:
        assembly=any_given("assembly")
    resources:
        parallel_downloads=1,
        genomepy_downloads=1,
    priority: 10
    retries: 2
    script:
        f"{config['rule_dir']}/../scripts/genomepy/get_genome.py"


rule get_genome_support_files:
    """
    Generate supporting files for a genome.
    """
    input:
        rules.get_genome.output,
    output:
        index=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.fai", **config),
        sizes=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
        gaps=expand("{genome_dir}/{{assembly}}/{{assembly}}.gaps.bed", **config),
    wildcard_constraints:
        assembly=".+",
    script:
        f"{config['rule_dir']}/../scripts/genomepy/get_genome_support.py"


rule get_genome_blacklist:
    """
    Download a genome blacklist for a genome, if it exists, through genomepy.
    """
    input:
        rules.get_genome.output,
        ancient(rules.get_genome_support_files.output)
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.blacklist.bed", **config),
    log:
        expand("{log_dir}/get_genome/{{assembly}}.blacklist.log", **config),
    wildcard_constraints:
        assembly=any_given("assembly")
    resources:
        parallel_downloads=1,
        genomepy_downloads=1,
    priority: 10
    retries: 2
    script:
        f"{config['rule_dir']}/../scripts/genomepy/get_genome_blacklist.py"


rule get_genome_annotation:
    """
    Download a gene annotation through genomepy.
    """
    input:
        rules.get_genome.output,
        ancient(rules.get_genome_support_files.output)
    output:
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        bed=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.bed", **config),
    log:
        expand("{log_dir}/get_genome/{{assembly}}.annotation.log", **config),
    params:
        providers=PROVIDERS,
        provider=config.get("provider"),
        genomes_dir=config["genome_dir"],
    wildcard_constraints:
        assembly=any_given("assembly")
    resources:
        parallel_downloads=1,
        genomepy_downloads=1,
    priority: 10
    retries: 2
    script:
        f"{config['rule_dir']}/../scripts/genomepy/get_genome_annotation.py"


rule extend_genome:
    """
    Append given file(s) to genome
    """
    input:
        genome=rules.get_genome.output,
        extension=config.get("custom_genome_extension", []),
    output:
        genome=expand(
            "{genome_dir}/{{assembly}}{custom_assembly_suffix}/{{assembly}}{custom_assembly_suffix}.fa", **config
        ),
    wildcard_constraints:
        assembly=any_given("assembly")
    message: EXPLAIN["custom_extension"]
    shell:
        """
        # extend the genome.fa
        cp {input.genome} {output.genome}

        for FILE in {input.extension}; do
            cat $FILE >> {output.genome}
        done
        """


rule extend_genome_blacklist:
    """
    Copy blacklist to the custom genome directory
    """
    input:
        rules.get_genome_blacklist.output,
    output:
        expand("{genome_dir}/{{assembly}}{custom_assembly_suffix}/{{assembly}}{custom_assembly_suffix}.blacklist.bed", **config),
    wildcard_constraints:
        assembly=any_given("assembly")
    shell:
        """
        cp {input} {output}
        """


rule extend_genome_annotation:
    """
    Append given file(s) to genome annotation
    """
    input:
        gtf=rules.get_genome_annotation.output.gtf,
        extension=config.get("custom_annotation_extension", []),
    output:
        gtf=expand("{genome_dir}/{{assembly}}{custom_assembly_suffix}/{{assembly}}{custom_assembly_suffix}.annotation.gtf", **config),
        bed=expand("{genome_dir}/{{assembly}}{custom_assembly_suffix}/{{assembly}}{custom_assembly_suffix}.annotation.bed", **config),
        gp=temp(
            expand("{genome_dir}/{{assembly}}{custom_assembly_suffix}/{{assembly}}{custom_assembly_suffix}.annotation.gp", **config)
        ),
    wildcard_constraints:
        assembly=any_given("assembly")
    message: EXPLAIN["custom_extension"]
    shell:
        """
        # extend the genome.annotation.gtf
        cp {input.gtf} {output.gtf}

        for FILE in {input.extension}; do
            cat $FILE >> {output.gtf}
        done

        # generate an extended genome.annotation.bed
        gtfToGenePred {output.gtf} {output.gp}
        genePredToBed {output.gp} {output.bed}
        """


rule gene_id2name:
    """
    Parse the gtf file to generate a gene_id to gene_name conversion table.
    """
    input:
        rules.get_genome_annotation.output.gtf,
    output:
        expand("{genome_dir}/{{assembly}}/gene_id2name.tsv", **config),
    script:
        f"{config['rule_dir']}/../scripts/gene_id2name.py"


rule get_effective_genome_size:
    """
    Get the effective genome size for a kmer length. Some tools (e.g. macs2) require
    an estimation of the effective genome size to better estimate how (un)likely it
    is to have a certain number of reads on a position. The actual genome size is
    not the best indication in these cases, since reads in repetitive regions
    (for a certain kmer length) can not possible align on some locations.
    """
    input:
        rules.get_genome.output,
    output:
        expand("{genome_dir}/{{assembly}}/genome_sizes/kmer_{{kmer_size}}.genome_size", **config),
    conda:
        "../envs/khmer.yaml"
    log:
        expand("{log_dir}/get_genome_size/{{assembly}}_{{kmer_size}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/get_genome_size/{{assembly}}_{{kmer_size}}.benchmark.txt", **config)[0]
    shell:
        """
        unique-kmers.py {input} -k {wildcards.kmer_size} --quiet 2>&1 | grep -P -o '(?<=\.fa: ).*' > {output} 2> {log}
        """
