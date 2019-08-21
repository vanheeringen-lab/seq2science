# the filetypes genomepy will download
config['genome_types'] = ['fa', 'fa.fai', "fa.sizes", "annotation.bed.gz", "annotation.gff.gz", "annotation.gtf"]

rule get_genome:
    """
    Download a genome through genomepy.

    Automatically turns on/off plugins.
    """
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
        dir=config['genome_dir'],
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf.gz", **config),
    conda:
        "../envs/get_genome.yaml"
    shell:
        """
        active_plugins=$(genomepy config show | grep -Po '(?<=- ).*' | paste -s -d, -)
        trap "genomepy plugin enable {{$active_plugins,}} > {log} 2>&1" EXIT

        genomepy plugin disable {{blacklist,bowtie2,bwa,gaps,gmap,hisat2,minimap2}} >> {log} 2>&1

        genomepy install --genome_dir {params.dir} {wildcards.assembly} UCSC --annotation    >> {log} 2>&1 ||
        genomepy install --genome_dir {params.dir} {wildcards.assembly} NCBI --annotation    >> {log} 2>&1 ||
        genomepy install --genome_dir {params.dir} {wildcards.assembly} Ensembl --annotation >> {log} 2>&1
        
        gunzip {params.gtf} >> {log} 2>&1
        """


rule get_transcripts:
    """
    Combine genome.fasta and annotation.gtf to create transcripts.fasta using gffread.

    Prepare genome.fasta if needed.
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
    priority: 1
    params:
        path=conda_path("../../envs/get_genome.yaml"),
        purgedfa=expand("{genome_dir}/{{assembly}}/{{assembly}}.purged.fa", **config)
    run:
        conda_gffread = os.path.join(params.path, "bin", "gffread")

        # Check if fasta has dirty formatting
        with open(input.fa[0], 'r') as fa:
            for line in fa:
                if line.startswith('>'):
                    line = line.strip().split(' ')
                    clean = True if len(line) == 1 else False
                    break

        if clean:
            shell(conda_gffread + " -w {output} -g {input.fa} {input.gtf} >> {log} 2>&1")

        # purge the dirty formatting from the fasta file, then use this in gffread
        else:
            # get location identifiers from the gtf (always the first element)
            locs = []
            with open(input.gtf[0], 'r') as gtf:
            #with open(params.gtf[0], 'r') as gtf:
                for line in gtf:
                    line = line.strip().split('\t')[0]
                    if line not in locs:
                        locs.append(line)

            # determine which element in the line contains the location identifier (checking the first only)
            with open(input.fa[0], 'r') as fa:
                for line in fa:
                    if line.startswith('>'):
                        line = line.strip(">\n").split(' ')
                        for loc in locs:
                            try:
                                element = line.index(loc)
                                break
                            except:
                                continue
                        break

            # extract the identifier from the header lines, keep the rest as-is.
            with open(input.fa[0], 'r') as fa:
                with open(params.purgedfa[0], 'w') as out:
                    for line in fa:
                        if line.startswith('>'):
                            identifier = line.strip(">").split(' ')[element]
                            line = '>' + identifier + '\n'
                        out.write(''.join(line))

            shell(conda_gffread + " -w {output} -g {params.purgedfa} {input.gtf} >> {log} 2>&1")
            