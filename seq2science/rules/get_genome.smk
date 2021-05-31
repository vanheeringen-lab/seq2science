localrules: extend_genome, get_genome_support_files, unzip_annotation

rule get_genome:
    """
    Download a genome through genomepy.
    """
    output:
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa", **config),
    log:
        expand("{log_dir}/get_genome/{{raw_assembly}}.genome.log", **config),
    benchmark:
        expand("{benchmark_dir}/get_genome/{{raw_assembly}}.genome.benchmark.txt", **config)[0]
    message: explain_rule("get_genome")
    params:
        providers=providers,
        genome_dir=config["genome_dir"]
    resources:
        parallel_downloads=1,
    priority: 1
    script:
        f"{config['rule_dir']}/../scripts/get_genome.py"


rule get_genome_blacklist:
    """
    Download a genome through genomepy.
    
    Also download a blacklist if it exists.
    """
    input:
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa", **config),
        # get_genome_support_files must have finished
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa.fai", **config),
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa.sizes", **config),
    output:
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.blacklist.bed", **config),
    log:
        expand("{log_dir}/get_genome/{{raw_assembly}}.blacklist.log", **config),
    params:
        genome_dir=config["genome_dir"]
    resources:
        parallel_downloads=1,
    priority: 1
    script:
        f"{config['rule_dir']}/../scripts/get_genome.py"


rule get_genome_annotation:
    """
    Download a gene annotation through genomepy.
    """
    input:
        # these should be ancient(), but snakemake simply ignores ancient items (bug)
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa", **config),
        # get_genome_support_files must have finished
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa.fai", **config),
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa.sizes", **config),
    output:
        gtf=expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.annotation.gtf.gz", **config),
        bed=expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.annotation.bed.gz", **config),
    log:
        expand("{log_dir}/get_annotation/{{raw_assembly}}.genome.log", **config),
    benchmark:
        expand("{benchmark_dir}/get_annotation/{{raw_assembly}}.genome.benchmark.txt", **config)[0]
    resources:
        parallel_downloads=1,
    params:
        providers=providers,
        genome_dir=config["genome_dir"]
    priority: 1
    script:
        f"{config['rule_dir']}/../scripts/get_genome_annotation.py"


rule extend_genome:
    """
    Append given file(s) to genome
    """
    input:
        genome=expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa", **config),
        extension=config.get("custom_genome_extension", []),
    output:
        genome=expand("{genome_dir}/{{raw_assembly}}{custom_assembly_suffix}/{{raw_assembly}}{custom_assembly_suffix}.fa", **config),
    message: explain_rule("custom_extension")
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
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.blacklist.bed", **config),
    output:
        expand("{genome_dir}/{{raw_assembly}}{custom_assembly_suffix}/{{raw_assembly}}{custom_assembly_suffix}.blacklist.bed", **config),
    shell:
        """cp {input} {output}"""


rule extend_genome_annotation:
    """
    Append given file(s) to genome annotation
    """
    input:
        gtf=expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.annotation.gtf", **config),
        extension=config.get("custom_annotation_extension", [])
    output:
        gtf=expand("{genome_dir}/{{raw_assembly}}{custom_assembly_suffix}/{{raw_assembly}}{custom_assembly_suffix}.annotation.gtf", **config),
        bed=expand("{genome_dir}/{{raw_assembly}}{custom_assembly_suffix}/{{raw_assembly}}{custom_assembly_suffix}.annotation.bed", **config),
        gp=temp(expand("{genome_dir}/{{raw_assembly}}{custom_assembly_suffix}/{{raw_assembly}}{custom_assembly_suffix}.annotation.gp", **config)),
    message: explain_rule("custom_extension")
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


rule get_genome_support_files:
    """
    Generate supporting files for a genome.
    """
    input:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.fai", **config),
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
        expand("{genome_dir}/{{assembly}}/{{assembly}}.gaps.bed", **config),
    params:
        genome_dir=config["genome_dir"]
    script:
        f"{config['rule_dir']}/../scripts/genome_support.py"


rule gene_id2name:
    """
    Parse the gtf file to generate a gene_id to gene_name conversion table.
    """
    input:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
    output:
        expand("{genome_dir}/{{assembly}}/gene_id2name.tsv", **config),
    run:
        def can_convert():
            """check if we can make a conversion table at all"""
            with open(input[0]) as gtf:
                for n, line in enumerate(gtf):
                    line = line.lower()
                    if "gene_id" in line and "gene_name" in line:
                        return True
                    if n > 100:
                        break
                return False

        if not can_convert():
            with open(output[0], "w") as out:
                out.write("assembly does not contain both gene_ids and gene_names\n")
        else:

            # loop over the gtf and store the conversion in the table
            table = dict()
            with open(input[0]) as gtf:
                for line in gtf:
                    try:
                        attributes = line.split("\t")[8].split(";")
                        id = name = None
                        for attribute in attributes:
                            attribute = attribute.strip()
                            if attribute.lower().startswith("gene_id"):
                                id = attribute.split(" ")[1].strip('"')
                            if attribute.lower().startswith("gene_name"):
                                name = attribute.split(" ")[1].strip('"')
                        if id and name:
                            table[id] = name
                    except IndexError:
                        # skip lines that are too short/misformatted
                        continue

            # save the dict
            with open(output[0], "w") as out:
                for k,v in table.items():
                    out.write(f"{k}\t{v}\n")


rule unzip_annotation:
    """
    Unzip (b)gzipped files.
    """
    input:
        "{filepath}.gz"
    output:
        "{filepath}"
    wildcard_constraints:
        filepath=".*(\.annotation)(\.gtf|\.bed)(?<!\.gz)$"  # filepath may not end with ".gz"
    priority: 1
    run:
        import genomepy.utils
        genomepy.utils.gunzip_and_name(input[0])
