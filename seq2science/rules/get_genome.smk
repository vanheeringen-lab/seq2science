import os
from functools import lru_cache
import contextlib

import genomepy


# the filetypes genomepy will download
config['genome_types'] = ['fa', 'fa.fai', "fa.sizes", "gaps.bed"]
config['genomepy_temp'] = ["annotation.bed.gz", "annotation.gff.gz"]

# add annotation to the expected output if it is required
if 'rna_seq' in get_workflow() or config['aligner'] == 'star':
    config['genome_types'].append("annotation.gtf")


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
        expand("{genome_dir}/{{assembly}}/{{assembly}}.{genome_types}", **config)
    log:
        expand("{log_dir}/get_genome/{{assembly}}.genome.log", **config)
    benchmark:
        expand("{benchmark_dir}/get_genome/{{assembly}}.genome.benchmark.txt", **config)[0]
    resources:
        parallel_downloads=1
    priority: 1
    params:
        dir=config['genome_dir'],
        provider=config.get('provider', None),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        temp=expand("{genome_dir}/{{assembly}}/{{assembly}}.{genomepy_temp}", **config)
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


rule get_transcripts:
    """
    Generate transcripts.fasta using gffread.
    
    Requires genome.fa and annotation.gtf (with matching chromosome/scaffold names)
    """
    input:
        fa=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config)
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config)
    log:
        expand("{log_dir}/get_genome/{{assembly}}.transcripts.log", **config)
    benchmark:
        expand("{benchmark_dir}/get_genome/{{assembly}}.transcripts.benchmark.txt", **config)[0]
    conda:
        "../envs/get_genome.yaml"
    priority: 1
    shell:
        "gffread -w {output} -g {input.fa} {input.gtf} >> {log} 2>&1"


rule decoy_transcripts:
    """
    Generate decoy_transcripts.txt for Salmon indexing  
    
    script source: https://github.com/COMBINE-lab/SalmonTools
    """
    input:
        genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
        transcripts=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
    output:
        expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config)
    params:
        script=f"{config['rule_dir']}/../scripts/generateDecoyTranscriptome.sh"
    log:
        expand("{log_dir}/get_genome/{{assembly}}.decoy_transcripts.log", **config)
    benchmark:
        expand("{benchmark_dir}/get_genome/{{assembly}}.decoy_transcripts.benchmark.txt", **config)[0]
    threads: 40
    resources:
        mem_gb=65
    conda:
        "../envs/decoy.yaml"
    priority: 1
    shell:
         ("cpulimit --include-children -l {threads}00 -- " if config.get("cpulimit", True) else " ") +
         "sh {params.script} -j {threads} -g {input.genome} -a {input.gtf} -t {input.transcripts} -o $(dirname {output}) > {log} 2>&1"


@lru_cache(maxsize=None)
def has_annotation(assembly):
    """
    returns True/False on whether or not the assembly has an annotation.
    """
    # check if genome providd by user or already downloaded, if so check if the annotation came along
    if all(os.path.exists(f"{config['genome_dir']}/{assembly}.{extension}") for extension in config["genome_types"]):
        return os.path.exists(f"{config['genome_dir']}/{assembly}.annotation.gtf")

    if "provider" in config:
        providers = config["provider"]
    else:
        providers = ["Ensembl", "UCSC", "NCBI"]

    # check if we expect an annotation
    # we do not want genomepy outputting to us that it is downloading stuff
    with open(os.devnull, 'w') as null:
        with contextlib.redirect_stdout(null), contextlib.redirect_stderr(null):
            for provider in providers:
                p = genomepy.ProviderBase.create(provider)
                if assembly in p.genomes and \
                    p.get_genome_download_link(assembly) is not None:
                    return p.get_annotation_download_link(assembly) is not None

    # no download link found for assembly
    return False
