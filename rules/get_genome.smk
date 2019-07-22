# the filetypes genomepy will download
config['genome_types'] = ['fa', 'fa.fai', 'fa.sizes', 'gaps.bed']

rule get_genome:
    # """
    # Download the genomes, specified in config:samples, through genomepy.
    # """
    output:
        expand("{genome_dir}/{{assemb}}/{{assemb}}.{genome_types}", **config)
    log:
        expand("logs/get_genome/{{assemb}}.log", **config)
    threads: 1
    resources:
        parallel_downloads=1
    priority: 1
    conda:
        "../envs/get_genome.yaml"
    shell:
        "genomepy install {wildcards.assemb} UCSC    >  {log} 2>&1 || "
        "genomepy install {wildcards.assemb} NCBI    >> {log} 2>&1 || "
        "genomepy install {wildcards.assemb} Ensembl >> {log} 2>&1 "
