import os
import time
import contextlib
import genomepy
from filelock import FileLock
from functools import lru_cache


# the filetypes genomepy will download
config["genome_types"] = ["fa", "fa.fai", "fa.sizes", "gaps.bed"]
config["genomepy_temp"] = ["annotation.bed.gz", "annotation.gff.gz"]

# add annotation to the expected output if it is required
if "rna_seq" in get_workflow() or config["aligner"] == "star":
    config["genome_types"].append("annotation.gtf")


# TODO: return to checkpoint get_genome when checkpoints are stable
#  1) does the trackhub input update? 2) does ruleorder work?
rule get_genome:
    """
    Download a genome through genomepy.
    Additionally downloads the gene annotation if required downstream.

    If assemblies with the same name can be downloaded from multiple providers, 
    a provider may be specified in the config (example: provider: NCBI). Otherwise,
    each provider will be tried in turn, stopping at the first success.

    Automatically turns on/off plugins.
    """
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.{genome_types}", **config),
    log:
        expand("{log_dir}/get_genome/{{assembly}}.genome.log", **config),
    benchmark:
        expand("{benchmark_dir}/get_genome/{{assembly}}.genome.benchmark.txt", **config)[0]
    message: explain_rule("get_genome")
    resources:
        parallel_downloads=1,
    priority: 1
    params:
        dir=config["genome_dir"],
        provider=config.get("provider", None),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        temp=expand("{genome_dir}/{{assembly}}/{{assembly}}.{genomepy_temp}", **config),
    conda:
        "../envs/get_genome.yaml"
    shell:
        """
        # turn off plugins and reset on exit. delete temp files on exit.
        active_plugins=$(genomepy config show | grep -Po '(?<=- ).*' | paste -s -d, -) || echo ""
        trap "genomepy plugin enable {{$active_plugins,}} >> {log} 2>&1; rm -f {params.temp}" EXIT
        genomepy plugin disable {{blacklist,bowtie2,bwa,star,gmap,hisat2,minimap2}} >> {log} 2>&1
        genomepy plugin enable {{blacklist,}} >> {log} 2>&1

        # download the genome and attempt to download the annotation and blacklist
        if [[ ! {params.provider} = None  ]]; then
            genomepy install --genomes_dir {params.dir} {wildcards.assembly} {params.provider} --annotation >> {log} 2>&1
        else
            genomepy install --genomes_dir {params.dir} {wildcards.assembly} Ensembl --annotation >> {log} 2>&1 ||
            genomepy install --genomes_dir {params.dir} {wildcards.assembly} UCSC    --annotation >> {log} 2>&1 ||
            genomepy install --genomes_dir {params.dir} {wildcards.assembly} NCBI    --annotation >> {log} 2>&1
        fi

        # unzip annotation if downloaded and gzipped
        if [ -f {params.gtf}.gz ]; then
            gunzip -f {params.gtf}.gz >> {log} 2>&1
        fi

        # if assembly has no annotation, or annotation has no genes, throw an warning
        if [ ! -f {params.gtf} ] || $(grep -q "No genes found" {log}); then
            echo '\nEmpty or no annotation found for {wildcards.assembly}.\n' | tee -a {log}

            # if an annotation is required, make it an error and exit.
            if $(echo {output} | grep -q annotation.gtf); then
                echo'\nSelect a different assembly or provide an annotation file manually.\n\n' | tee -a {log}
                exit 1
            fi
        fi
        """


@lru_cache(maxsize=None)
def has_annotation(assembly):
    """
    returns True/False on whether or not the assembly has an annotation.
    """
    # check if genome is provided by user or already downloaded, if so check if the annotation came along
    if all(os.path.exists(f"{config['genome_dir']}/{assembly}.{extension}") for extension in config["genome_types"]):
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
