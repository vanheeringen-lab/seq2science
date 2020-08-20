import contextlib
import genomepy


rule get_genome:
    """
    Download a genome through genomepy.
    
    Also download a blacklist if it exists.
    """
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
    log:
        expand("{log_dir}/get_genome/{{assembly}}.genome.log", **config),
    benchmark:
        expand("{benchmark_dir}/get_genome/{{assembly}}.genome.benchmark.txt", **config)[0]
    message: explain_rule("get_genome")
    resources:
        parallel_downloads=1,
    priority: 1
    run:
        with open(log[0], "w") as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):

                # select a provider with the annotation if possible
                a = providers[wildcards.assembly]["annotation"]
                g = providers[wildcards.assembly]["genome"]
                provider = g if a is None else a

                p = genomepy.ProviderBase.create(provider)
                p.download_genome(wildcards.assembly, config["genome_dir"])

                # try to download the blacklist
                genome = genomepy.Genome(wildcards.assembly, config["genome_dir"])
                plugins = genomepy.plugin.init_plugins()
                plugins["blacklist"].after_genome_download(genome)


rule get_genome_annotation:
    """
    Download a gene annotation through genomepy.
    """
    input:
        rules.get_genome.output,
    output:
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf.gz", **config),
        bed=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.bed.gz", **config),
    log:
        expand("{log_dir}/get_annotation/{{assembly}}.genome.log", **config),
    benchmark:
        expand("{benchmark_dir}/get_annotation/{{assembly}}.genome.benchmark.txt", **config)[0]
    resources:
        parallel_downloads=1,
    priority: 1
    run:
        with open(log[0], "w") as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):

                provider = providers[wildcards.assembly]["annotation"]
                p = genomepy.ProviderBase.create(provider)

                kwargs = dict()
                if provider == "ucsc" and p.get_annotation_download_link(
                        name=wildcards.assembly, kwargs={"ucsc_annotation_type": "ensembl"}):
                    kwargs = {"ucsc_annotation_type": "ensembl"}

                p.download_annotation(
                    name=wildcards.assembly,
                    genomes_dir=config["genome_dir"],
                    localname=wildcards.assembly,
                    kwargs=kwargs)

                # sanitize the annotations
                genome = genomepy.Genome(wildcards.assembly, config["genome_dir"])
                genomepy.utils.sanitize_annotation(genome)


rule get_genome_support_files:
    """
    Generate supporting files for a genome.
    """
    input:
        expand("{genome_dir}/{{assembly_}}/{{assembly_}}.fa", **config),
    output:
        expand("{genome_dir}/{{assembly_}}/{{assembly_}}.fa.fai", **config),
        expand("{genome_dir}/{{assembly_}}/{{assembly_}}.fa.sizes", **config),
        expand("{genome_dir}/{{assembly_}}/{{assembly_}}.gaps.bed", **config),
    run:
        genomepy.Genome(wildcards.assembly_, genomes_dir=config["genome_dir"])


rule add_spike_ins:
    """
    Create a "new" genome which includes the spike-in data.
    """
    input:
        genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        spike_in_fa=config.get("spike_in_fa", ""),
        spike_in_gtf=config.get("spike_in_gtf", "")
    output:
        genome=expand("{genome_dir}/{{assembly}}_SI/{{assembly}}_SI.fa", **config),
        gtf=expand("{genome_dir}/{{assembly}}_SI/{{assembly}}_SI.annotation.gtf", **config),
        bed=expand("{genome_dir}/{{assembly}}_SI/{{assembly}}_SI.annotation.bed", **config),
        gp=temp(expand("{genome_dir}/{{assembly}}_SI/{{assembly}}_SI.annotation.gp", **config)),
    message: explain_rule("add_spike_ins")
    shell:
        """
        # add spike ins to the genome.fa
        cp {input.genome} {output.genome}
        cat {input.spike_in_fa} >> {output.genome}

        # add spike ins to the genome.annotation.gtf
        cp {input.gtf} {output.gtf}
        cat {input.spike_in_gtf} >> {output.gtf}

        # generate spike-in genome.annotation.bed
        gtfToGenePred {output.gtf} {output.gp}
        genePredToBed {output.gp} {output.bed}
        """


# NOTE: if the workflow fails it tends to blame this rule.
# Set "debug: True" in the config, or comment out the rule to see the root cause.
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
