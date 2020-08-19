# import os
# import time
import contextlib
import genomepy
# from filelock import FileLock
# from functools import lru_cache





# @lru_cache(maxsize=None)
# def has_annotation(assembly):
#     """
#     Returns True/False on whether or not the assembly has an annotation.
#     """
#     # check if the annotation exists locally
#     gtf = os.path.join(config['genome_dir'], assembly, f"{assembly}.annotation.gtf")
#     if any(os.path.exists(file) for file in [gtf, f"{gtf}.gz"]):
#         return True
#
#     # check if the annotation can be downloaded
#     if provider_with_file("annotation", assembly):
#         return True
#
#     # # warn if an annotation is not needed, exit if it is.
#     # logger.info(
#     #     f"No annotation for assembly {assembly} can be downloaded. Another provider (and "
#     #     f"thus another assembly name) might have gene annotations.\n"
#     #     f"Find alternative assemblies with `genomepy search {assembly}`"
#     # )
#     # annotion_required = "rna_seq" in get_workflow() or config["aligner"] == "star"
#     # if annotion_required:
#     #     exit(1)
#     # time.sleep(3)
#     return False


# rule get_genome_provider:
#     """
#     Determine which provider to download the assembly from.
#
#     Uses the provider used before, or else the provider specified in the config (example: provider: NCBI).
#     If neither is found, each provider will be tried.
#     """
#     output:
#         temp(expand("{genome_dir}/{{assembly}}/provider.txt", **config)),
#     resources:
#         parallel_downloads=1,
#     priority: 2
#     run:
#         # genome_provider = provider_with_file("genome", wildcards.assembly)
#         # if not genome_provider:
#         #     logger.info(
#         #         f"Could not download assembly {wildcards.assembly}.\n"
#         #         f"Find alternative assemblies with `genomepy search {wildcards.assembly}`"
#         #     )
#
#         # annotation_required = "rna_seq" in get_workflow() or config["aligner"] == "star"
#         # annotation_provider = provider_with_file("annotation", wildcards.assembly)
#         # if annotation_provider:
#         #     provider = annotation_provider
#         # elif has_annotation(wildcards.assembly) or not annotation_required:
#         #     # assumes that if a local annotation is present that
#         #     # the user made sure it is in the same style
#         #     provider = provider_with_file("genome", wildcards.assembly)
#         # else:
#         #     logger.info(
#         #         f"No annotation for assembly {wildcards.assembly} can be downloaded. Another provider (and "
#         #         f"thus another assembly name) might have gene annotations.\n"
#         #         f"Find alternative assemblies with `genomepy search {wildcards.assembly}`"
#         #     )
#         #     exit(1)
#         annotation_provider = provider_with_file("annotation", wildcards.assembly)
#         if annotation_provider:
#             provider = annotation_provider
#         else:
#             provider = provider_with_file("genome", wildcards.assembly)
#
#         with open(output[0], "w") as out:
#             out.write(provider)


rule get_genome:
    """
    Download a genome through genomepy.
    
    Also download a blacklist if it exists.
    """
    # input:
    #     rules.get_genome_provider.output
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
                # provider = open(input[0]).read()
                provider = providers[wildcards.assembly]
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
        # rules.get_genome_provider.output,
        rules.get_genome.output,  # required for sanitizing
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
                # provider = open(input[0]).read()
                provider = providers[wildcards.assembly]
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

                # try to sanitize the annotations
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
    run:
        genomepy.utils.gunzip_and_name(input[0])
