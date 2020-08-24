import contextlib
import genomepy


# the filetypes genomepy will download
config["genome_types"] = ["fa", "fa.fai", "fa.sizes", "gaps.bed"]
config["genomepy_temp"] = ["annotation.gff.gz"]

# add annotation to the expected output if it is required
if "rna_seq" == get_workflow() or config["aligner"] == "star" or \
        "scrna_seq" == get_workflow():
    config["genome_types"].extend(["annotation.gtf", "annotation.bed"])


# TODO: return to checkpoint get_genome when checkpoints are stable
#  1) does the trackhub input update? 2) does ruleorder work?
rule get_genome:
    """
    Download a genome through genomepy.
    
    Also download a blacklist if it exists.
    """
    output:
        expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa", **config),
    log:
        expand("{log_dir}/get_genome/{{raw_assembly}}.genome.log", **config),
    benchmark:
        expand("{benchmark_dir}/get_genome/{{raw_assembly}}.genome.benchmark.txt", **config)[0]
    message: explain_rule("get_genome")
    resources:
        parallel_downloads=1,
    priority: 1
    run:
        with open(log[0], "w") as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):

                # select a provider with the annotation if possible
                a = providers[wildcards.raw_assembly]["annotation"]
                g = providers[wildcards.raw_assembly]["genome"]
                provider = g if a is None else a

                p = genomepy.ProviderBase.create(provider)
                p.download_genome(wildcards.raw_assembly, config["genome_dir"])

                # try to download the blacklist
                genome = genomepy.Genome(wildcards.raw_assembly, config["genome_dir"])
                plugins = genomepy.plugin.init_plugins()
                plugins["blacklist"].after_genome_download(genome)


rule get_genome_annotation:
    """
    Download a gene annotation through genomepy.
    """
    input:
        rules.get_genome.output,
    output:
        gtf=expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.annotation.gtf.gz", **config),
        bed=expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.annotation.bed.gz", **config),
    log:
        expand("{log_dir}/get_annotation/{{raw_assembly}}.genome.log", **config),
    benchmark:
        expand("{benchmark_dir}/get_annotation/{{raw_assembly}}.genome.benchmark.txt", **config)[0]
    resources:
        parallel_downloads=1,
    priority: 1
    run:
        with open(log[0], "w") as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):

                provider = providers[wildcards.raw_assembly]["annotation"]
                p = genomepy.ProviderBase.create(provider)

                kwargs = dict()
                if provider == "ucsc" and p.get_annotation_download_link(
                        name=wildcards.raw_assembly, kwargs={"ucsc_annotation_type": "ensembl"}):
                    kwargs = {"ucsc_annotation_type": "ensembl"}

                p.download_annotation(
                    name=wildcards.raw_assembly,
                    genomes_dir=config["genome_dir"],
                    localname=wildcards.raw_assembly,
                    kwargs=kwargs)

                # sanitize the annotations
                genome = genomepy.Genome(wildcards.raw_assembly, config["genome_dir"])
                genomepy.utils.sanitize_annotation(genome)


rule extend_genome:
    """
    Append given file(s) to genome
    """
    input:
        genome=expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.fa", **config),
        extension=config.get("custom_genome_extension", []),
    output:
        genome=expand("{genome_dir}/{{raw_assembly}}_custom/{{raw_assembly}}_custom.fa", **config),
    message: explain_rule("custom_extension")
    shell:
        """
        # extend the genome.fa
        cp {input.genome} {output.genome}
        
        for FILE in {input.extension}; do
            cat $FILE >> {output.genome}
        done
        """


rule extend_genome_annotation:
    """
    Append given file(s) to genome annotation
    """
    input:
        gtf=expand("{genome_dir}/{{raw_assembly}}/{{raw_assembly}}.annotation.gtf", **config),
        extension=config.get("custom_annotation_extension", [])
    output:
        gtf=expand("{genome_dir}/{{raw_assembly}}_custom/{{raw_assembly}}_custom.annotation.gtf", **config),
        bed=expand("{genome_dir}/{{raw_assembly}}_custom/{{raw_assembly}}_custom.annotation.bed", **config),
        gp=temp(expand("{genome_dir}/{{raw_assembly}}_custom/{{raw_assembly}}_custom.annotation.gp", **config)),
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
    run:
        genomepy.Genome(wildcards.assembly, genomes_dir=config["genome_dir"])


# NOTE: if the workflow fails it tends to blame this rule.
# Set "debug: True" in the config to see the root cause.
if not config.get("debug"):
    rule unzip_file:
        """
        Unzip (b)gzipped files.
        """
        input:
            "{filepath}.gz"
        output:
            "{filepath}"
        wildcard_constraints:
            filepath=".*(?<!\.gz)$"  # filepath may not end with ".gz"
        priority: 1
        run:
            genomepy.utils.gunzip_and_name(input[0])
