# the filetypes genomepy will download
config['genome_types'] = ['fa', 'fa.fai', "fa.sizes", "gaps.bed"]
config['genomepy_temp'] = ["annotation.bed.gz", "annotation.gff.gz"]
# add annotation to the expected output if it is required
if 'rna_seq' in get_workflow() or config['aligner'] == 'star':
    config['genome_types'].append("annotation.gtf")

checkpoint get_genome:
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
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf.gz", **config),
        temp=expand("{genome_dir}/{{assembly}}/{{assembly}}.{genomepy_temp}", **config)
    conda:
        "../envs/get_genome.yaml"
    shell:
        """
        # turn off plugins and reset on exit. delete temp files on exit.
        active_plugins=$(genomepy config show | grep -Po '(?<=- ).*' | paste -s -d, -)
        trap "genomepy plugin enable {{$active_plugins,}} >> {log} 2>&1; rm -f {params.temp}" EXIT
        genomepy plugin disable {{blacklist,bowtie2,bwa,star,gmap,hisat2,minimap2}} >> {log} 2>&1
        
        # download the genome and attempt to download the annotation
        if [[ ! {params.provider} = None  ]]; then
            genomepy install --genome_dir {params.dir} {wildcards.assembly} {params.provider} --annotation >> {log} 2>&1
        else
            genomepy install --genome_dir {params.dir} {wildcards.assembly} Ensembl --annotation >> {log} 2>&1 ||
            genomepy install --genome_dir {params.dir} {wildcards.assembly} UCSC    --annotation >> {log} 2>&1 ||
            genomepy install --genome_dir {params.dir} {wildcards.assembly} NCBI    --annotation >> {log} 2>&1
        fi
        
        # unzip annotation if downloaded, warn if required but empty
        if [ -f {params.gtf} ]; then
            gunzip {params.gtf} >> {log} 2>&1
        elif echo {output} | grep -q annotation.gtf; then
            echo '\nAnnotation for {wildcards.assembly} contains no genes. Select a different assembly or provide an annotation file manually.\n\n' > {log}
            exit 1
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

        with open(input.sizefile[0], 'r') as sizes:
            for line in sizes:
                fa_id = line.split('\t')[0]
                if fa_id == gtf_id:
                    shell("echo 'Genome and annotation have matching chromosome/scaffold names!' >> {log}")
                    shell('ln {input.gtf} {output}')
                    break
            else:
                # generate a gtf with matching scaffold/chromosome IDs
                shell("echo 'Genome and annotation do not have matching chromosome/scaffold names! Creating matching gtf...' >> {log}")

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
                with open(input.gtf[0], 'r') as oldgtf, open(output[0], 'w') as newgtf:
                    for line in oldgtf:
                        line = line.split('\t')
                        line[0] = ids[line[0]]
                        line = '\t'.join(line)
                        newgtf.write(line)

                shell("echo '\nCorrected GTF creation completed.' >> {log}")


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
    conda:
        "../envs/get_genome.yaml"
    priority: 1
    shell:
        "gffread -w {output} -g {input.fa} {input.gtf} >> {log} 2>&1"


rule decoy_transcripts:
    """
    Generate decoy_transcripts.txt
    
    script source: https://github.com/COMBINE-lab/SalmonTools
    """
    input:
        script="../../scripts/generateDecoyTranscriptome.sh",
        genome=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
        gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.gtf", **config),
        transcripts=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
    output:
        expand("{genome_dir}/{{assembly}}/decoy_transcripts/decoys.txt", **config)
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
         """
         outdir=$(dirname {output})
         
         cpulimit --include-children -l {threads}00 --\
         sh {input.script} -j {threads} -g {input.genome} -a {input.gtf} -t {input.transcripts} -o $outdir > {log} 2>&1
         """
