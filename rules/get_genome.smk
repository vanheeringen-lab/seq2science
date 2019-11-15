# the filetypes genomepy will download
config['genome_types'] = ['fa', 'fa.fai', "fa.sizes"]
# intermediate filetypes genomepy may download
config['genomepy_temp'] = []

# transcript & annotation related filetypes
if config['aligner'] in ['salmon', 'star']:
    config['genome_types'].append("annotation.gtf")
    config['genomepy_temp'].extend(["annotation.bed.gz", "annotation.gff.gz"])


checkpoint get_genome:
    """
    Download a genome through genomepy.
    
    Additionally downloads the gene annotation if required downstream.

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
        annotation=lambda wildcards, output: '--annotation' if any('annotation.gtf' in file for file in output) else '',
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf.gz", **config),
        temp=expand("{genome_dir}/{{assembly}}/{{assembly}}.{genomepy_temp}", **config)
    conda:
        "../envs/get_genome.yaml"
    shell:
        """
        active_plugins=$(genomepy config show | grep -Po '(?<=- ).*' | paste -s -d, -)
        trap "genomepy plugin enable {{$active_plugins,}} >> {log} 2>&1; rm -f {params.temp}" EXIT

        genomepy plugin disable {{bowtie2,bwa,gaps,gmap,hisat2,minimap2}} >> {log} 2>&1
        
        genomepy install --genome_dir {params.dir} {wildcards.assembly} UCSC    {params.annotation} >> {log} 2>&1 ||
        genomepy install --genome_dir {params.dir} {wildcards.assembly} NCBI    {params.annotation} >> {log} 2>&1 ||
        genomepy install --genome_dir {params.dir} {wildcards.assembly} Ensembl {params.annotation} >> {log} 2>&1
        
        if [ {params.annotation} ]; then
            if [ -f {params.gtf} ]; then
                gunzip {params.gtf} >> {log} 2>&1
            else
                echo 'Annotation for {wildcards.assembly} contains no genes. Select a different assembly or provide an annotation file manually.\n\n' > {log}
                exit 1
            fi
        fi
        """


rule get_annotation:
    """
    Matches the chromosome/scaffold names in annotation.gtf to those in the genome.fa.
    """
    input:
        fa=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
        sizefile= expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config)
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.gtf", **config)
    log:
        expand("{log_dir}/get_genome/{{assembly}}.annotation.log", **config)
    benchmark:
        expand("{benchmark_dir}/get_genome/{{assembly}}.annotation.benchmark.txt", **config)[0]
    priority: 1
    run:
        # Check if genome and annotation have matching chromosome/scaffold names
        with open(input.gtf[0], 'r') as gtf:
            for line in gtf:
                if not line.startswith('#'):
                    gtf_id = line.split('\t')[0]
                    break

        matching = False
        with open(input.sizefile[0], 'r') as sizes:
            for line in sizes:
                fa_id = line.split('\t')[0]
                if fa_id == gtf_id:
                    matching = True
                    break

        if matching:
            shell("echo 'genome and annotation have matching chromosome/scaffold names!' >> {log}")
            shell('mv {input.gtf} {output}')

        else:
            # generate a gtf with matching scaffold/chromosome IDs
            shell("echo 'genome and annotation do not match, renaming gtf...' >> {log}")

            # determine which element in the genome.fasta's header contains the location identifiers used in the annotation.gtf
            header = []
            with open(input.fa[0], 'r') as fa:
                for line in fa:
                    if line.startswith('>'):
                        header = line.strip(">\n").split(' ')
                        break

            with open(input.gtf[0], 'r') as gtf:
                for line in gtf:
                    if not line.startswith('#'):
                        loc_id = line.strip().split('\t')[0]
                        try:
                            element = header.index(loc_id)
                            break
                        except:
                            continue

            # build a conversion table
            ids = {}
            with open(input.fa[0], 'r') as fa:
                for line in fa:
                    if line.startswith('>'):
                        line = line.strip(">\n").split(' ')
                        if line[element] not in ids.keys():
                            ids.update({line[element] : line[0]})

            # rename the location identifier in the gtf (using the conversion table)
            with open(input.gtf[0], 'r') as oldgtf:
                with open(output[0], 'w') as newgtf:
                    for line in oldgtf:
                        line = line.split('\t')
                        line[0] = ids[line[0]]
                        line = '\t'.join(line)
                        newgtf.write(line)

            shell("echo '\nRenaming successful! Deleting mismatched gtf file.' >> {log}")
            shell('rm -f {input.gtf}')


rule get_transcripts:
    """
    Generate transcripts.fasta using gffread.
    
    Requires genome.fa and annotation.gtf (with matching chromosome/scaffold names)
    """
    input:
        fa=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.gtf", **config)
    output:
        expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config)
    log:
        expand("{log_dir}/get_genome/{{assembly}}.transcripts.log", **config)
    benchmark:
        expand("{benchmark_dir}/get_genome/{{assembly}}.transcripts.benchmark.txt", **config)[0]
    priority: 1
    conda:
        "../envs/get_genome.yaml"
    shell:
        "gffread -w {output} -g {input.fa} {input.gtf} >> {log} 2>&1"
