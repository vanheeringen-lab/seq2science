import os
import time
import contextlib
import genomepy
from filelock import FileLock
from functools import lru_cache


rule get_genome_provider:
    """
    Determine which provider to download the assembly from.
    
    Uses the provider used before, or else the provider specified in the config (example: provider: NCBI).
    If neither is found, each provider will be tried.
    """
    output:
        temp(expand("{genome_dir}/{{assembly}}/provider.txt", **config)),
    resources:
        parallel_downloads=1,
    priority: 2
    run:
        readme_file = os.path.join(os.path.dirname(output[0]), "README.txt")
        readme_provider = None
        if os.path.exists(readme_file):
            metadata, _ = genomepy.utils.read_readme(readme_file)
            readme_provider = metadata.get("provider", "").lower()

        if readme_provider in ["ensembl", "ucsc", "ncbi"]:
            providers = [readme_provider]
        elif config.get("provider"):
            providers = [config["provider"].lower()]
        else:
            providers = list(genomepy.functions.list_available_providers())
            providers = providers[:-1]  # remove "url"

        annotion_required = "rna_seq" in get_workflow() or config["aligner"] == "star"
        for provider in providers:
            p = genomepy.ProviderBase.create(provider)
            if wildcards.assembly in p.genomes:
                annotation_found = p.get_annotation_download_link(wildcards.assembly)
                genome_found = p.get_genome_download_link(wildcards.assembly)
                if genome_found and annotation_found:
                    break
                elif genome_found and not annotion_required:
                    print("warning")
                    break
        else:
            print("error")
            exit(1)

        with open(output[0], "w") as out:
            out.write(provider)


rule get_genome:
    """
    Download a genome through genomepy.
    
    Also download a blacklist if it exists
    """
    input:
        rules.get_genome_provider.output
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
        provider = open(input[0]).read()
        p = genomepy.ProviderBase.create(provider)
        p.download_genome(wildcards.assembly, config["genome_dir"])

        # try to download the blacklist
        genome = genomepy.Genome(wildcards.assembly, config["genome_dir"])
        plugins = genomepy.plugin.init_plugins()
        plugins["blacklist"].after_genome_download(genome)


rule get_annotation:
    """
    Download a gene annotation through genomepy.
    """
    input:
        rules.get_genome_provider.output,
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
        provider = open(input[0]).read()
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
        rules.get_genome.output
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.fai", **config),
        expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
        expand("{genome_dir}/{{assembly}}/{{assembly}}.gaps.bed", **config),
    run:
        genomepy.Genome(wildcards.assembly, genomes_dir=config["genome_dir"])


rule add_spike_ins:
    """
    
    """
    input:
        genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        spike_ins=config.get("spike_ins")
    output:
        genome=expand("{genome_dir}/{{assembly}}_SI/{{assembly}}_SI.fa", **config),
        gtf=expand("{genome_dir}/{{assembly}}_SI/{{assembly}}_SI.annotation.gtf", **config),
        bed=expand("{genome_dir}/{{assembly}}_SI/{{assembly}}_SI.annotation.bed", **config),
        gp=temp(expand("{genome_dir}/{{assembly}}_SI/{{assembly}}_SI.annotation.gp", **config)),
    run:
        # add spike ins to the genome.fa
        # TODO check
        with open(input.genome, "r").read() as infile, open(input.spike_ins, "r").read() as spike_in_file, open(output.genome, "w") as outfile:
            outfile.write(infile)
            outfile.write(spike_in_file)

        # TODO: add spike in gtf

        # generate a new bed file from gtf
        shell("gtfToGenePred {output.gtf} {output.gp}")
        shell("genePredToBed {output.gp} {output.bed}")


rule unzip_file:
    """
    unzip gzipped or bgzipped files
    
    (doesnt crash if the file was unzipped before executing)
    """
    input:
        expand("{{filepath}}.gz")
    output:
        expand("{{filepath}}")
    run:
        genomepy.utils.gunzip_and_name(input[0])


@lru_cache(maxsize=None)
def has_annotation(assembly):
    """
    returns True/False on whether or not the assembly has an annotation.
    """
    # check if genome is provided by user or already downloaded, if so check if the annotation came along
    if all(os.path.exists(f"{config['genome_dir']}/{assembly}/{assembly}.{extension}") for extension in config["genome_types"]):
        return os.path.exists(f"{config['genome_dir']}/{assembly}.annotation.gtf")

    if "provider" in config:
        providers = [config["provider"]]
    else:
        providers = ["Ensembl", "UCSC", "NCBI"]

    # check if we expect an annotation
    # we do not want genomepy outputting to us that it is downloading stuff
    with open(os.devnull, "w") as null:
        with contextlib.redirect_stdout(null), contextlib.redirect_stderr(null):
            for provider in providers:
                annotation_lock = os.path.expanduser(f'~/.config/seq2science/genomepy_{provider}_annotations.lock')
                prep_filelock(annotation_lock, 20)

                with FileLock(annotation_lock):
                    p = genomepy.ProviderBase.create(provider)
                    if assembly in p.genomes:
                        if p.get_annotation_download_link(assembly) is None:
                            logger.info(
                                f"No annotation for assembly {assembly} can be downloaded. Another provider (and "
                                f"thus another assembly name) might have gene annotations.\n"
                                f"Find alternative assemblies with `genomepy search {assembly}`"
                            )
                            time.sleep(2)
                            return False
                        else:
                            return True

    # no download link found for assembly
    return False
