# the filetypes genomepy will download
config['genome_types'] = ['fa', 'fa.fai', 'fa.sizes', 'gaps.bed']

rule get_genome:
    # """
    # Download the genomes, specified in config:samples, through genomepy.
    # """
    output:
        expand("{genome_dir}/{{assemb}}/{{assemb}}.{genome_types}", **config)
    log:
        expand("{log_dir}/get_genome/{{assemb}}.log", **config)
    threads: 1
    resources:
        parallel_downloads=1
    priority: 1
    conda:
        "../envs/get_genome.yaml"
    shell:
        f"""
        genomepy install --genome_dir {config['genome_dir']} {{wildcards.assemb}} UCSC    >  {{log}} 2>&1 ||
        genomepy install --genome_dir {config['genome_dir']} {{wildcards.assemb}} NCBI    >> {{log}} 2>&1 ||
        genomepy install --genome_dir {config['genome_dir']} {{wildcards.assemb}} Ensembl >> {{log}} 2>&1
        """
