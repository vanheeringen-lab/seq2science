# the filetypes genomepy will download
config['genome_types'] = ['fa', 'fa.fai', "fa.sizes"]

rule get_genome:
    # """
    # Download genomes
    # """
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.{genome_types}", **config)
    log:
        expand("{log_dir}/get_genome/{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/get_genome/{{assembly}}.benchmark.txt", **config)[0]
    resources:
        parallel_downloads=1
    priority: 1
    params:
        config['genome_dir']
    conda:
        "../envs/get_genome.yaml"
    shell:
        """
        active_plugins=$(genomepy config show | grep -Po '(?<=- ).*' | paste -s -d, -)
        trap "genomepy plugin enable {{$active_plugins}} > {log} 2>&1" 0

        genomepy plugin disable {{blacklist,bowtie2,bwa,gaps,gmap,hisat2,minimap2}} >> {log} 2>&1

        genomepy install --genome_dir {params} {wildcards.assembly} UCSC    >> {log} 2>&1 ||
        genomepy install --genome_dir {params} {wildcards.assembly} NCBI    >> {log} 2>&1 ||
        genomepy install --genome_dir {params} {wildcards.assembly} Ensembl >> {log} 2>&1
        """
