# the filetypes genomepy will download
config['genome_types'] = ['fa', 'fa.fai', "fa.sizes"]
if config['aligner'] == 'salmon': config['genome_types'].append("annotation.gtf")

# some of these intermediate files may be downloaded by genomepy
config['variable_temp_files'] = ["annotation.bed.gz", "annotation.gff.gz"]

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
        #annotation=lambda wildcards, output: 'annotation.gtf' in output[3], # '--annotation' if 'annotation.gtf' in output[3] else '',
        #gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf.gz", **config),
        gtf=lambda wildcards, output: expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf.gz", **config) if 'annotation.gtf' in output[3] else '',
        temp=expand("{genome_dir}/{{assembly}}/{{assembly}}.{variable_temp_files}", **config)
    conda:
        "../envs/get_genome.yaml"
    shell:
        """
        active_plugins=$(genomepy config show | grep -Po '(?<=- ).*' | paste -s -d, -)
        trap "genomepy plugin enable {{$active_plugins,}} >> {log} 2>&1; rm -f {params.temp}" EXIT

        genomepy plugin disable {{blacklist,bowtie2,bwa,gaps,gmap,hisat2,minimap2}} >> {log} 2>&1
        
        if [ {params.gtf} ]; then
            (genomepy install --genome_dir {params.dir} {wildcards.assembly} UCSC    --annotation >> {log} 2>&1 && [ -f {params.gtf} ]) ||
            (genomepy install --genome_dir {params.dir} {wildcards.assembly} NCBI    --annotation >> {log} 2>&1 && [ -f {params.gtf} ]) ||
            (genomepy install --genome_dir {params.dir} {wildcards.assembly} Ensembl --annotation >> {log} 2>&1 && [ -f {params.gtf} ]) ||
            (if cat {log} | grep -q "WARNING"; then echo '\n\n \
            Annotation for {wildcards.assembly} contains no genes. Select a different assembly or provide an annotation file manually.\
            \n\n' | tee -a {log}; fi && false)
            
            gunzip {params.gtf} >> {log} 2>&1
        else
            genomepy install --genome_dir {params.dir} {wildcards.assembly} UCSC    >> {log} 2>&1 ||
            genomepy install --genome_dir {params.dir} {wildcards.assembly} NCBI    >> {log} 2>&1 ||
            genomepy install --genome_dir {params.dir} {wildcards.assembly} Ensembl >> {log} 2>&1
        fi
        """

        # crashes the rul as desired
        # (genomepy install --genome_dir {params.dir} {wildcards.assembly} UCSC    --annotation  >> {log} 2>&1 && [ -f {params.gtf} ]) ||
        # (genomepy install --genome_dir {params.dir} {wildcards.assembly} NCBI    --annotation  >> {log} 2>&1 && [ -f {params.gtf} ]) ||
        # (genomepy install --genome_dir {params.dir} {wildcards.assembly} Ensembl --annotation  >> {log} 2>&1 && [ -f {params.gtf} ])

        # works but is ugly:
        # dir=$dirname {output[0]}
        # rm -rf $dir && mkdir $dir

        # tries to find the files before they are generated:
        # ls -d $(dirname {output[0]})/* | grep -v "$(cat {output})" | xargs rm -v

        # ls -d $(dirname $output)/* | grep -v $(basename $output) | xargs rm
        # ls -d $(dirname {output}[0])/* | grep -v {output} | xargs rm
        # ls -d $(dirname {output}[0])/* | grep -v $(basename {output}) | xargs rm

        #echo $(dirname {output}[0])
        #echo $(dirname /home/siebrenf/git/snakemake-workflows/workflows/alignment/xenTro9/xenTro9.fa)
        # rm -v $(dirname {output}[0]) !({output})
        # rm -rf $(dirname {output}[0])
        # rm -rf $(pwd)/.snakemake/shadow #deletes the whole shadow directory, which should always work, as a new dir is created upon starting a shadow rule
        # rm -rf $(ls -d $(pwd)/.snakemake/shadow/*) #works, but breaks the pipeline if the rule completes!


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

        else:
            # purge the dirty formatting from the fasta file, then use this in gffread
            # get location identifiers from the gtf (always the first element)
            locs = []
            with open(input.gtf[0], 'r') as gtf:
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
