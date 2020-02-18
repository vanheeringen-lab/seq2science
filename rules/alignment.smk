def get_reads(wildcards):
    if config.get('combine_replicates', '') == 'merge' and 'condition' in samples:
        if config['layout'].get(wildcards.sample, False) == "SINGLE":
            return expand("{trimmed_dir}/merged/{{sample}}_trimmed.{fqsuffix}.gz", **config)
        return sorted(expand("{trimmed_dir}/merged/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))
    else:
        if config['layout'].get(wildcards.sample, False) == "SINGLE":
            return expand("{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config)
        return sorted(expand("{trimmed_dir}/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))

def get_alignment_pipes():
    pipes = {pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])}
    if config.get('peak_caller', False) and 'genrich' in config['peak_caller']:
        pipes.add(pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.sambamba-queryname.pipe", **config)[0]))
    else:
        pipes.add(pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.{bam_sorter}-{bam_sort_order}.pipe", **config)[0]))

    return pipes


if config['aligner'] == 'bowtie2':
    rule bowtie2_index:
        """
        Make a genome index for bowtie2.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/bowtie2/", **config))
        log:
            expand("{log_dir}/bowtie2_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/bowtie2_index/{{assembly}}.benchmark.txt", **config)[0]
        priority: 1
        threads: 4
        conda:
            "../envs/bowtie2.yaml"
        params: config['index']
        shell:
            "bowtie2-build {params} --threads {threads} {input} {output}/{wildcards.assembly} > {log} 2>&1"


    rule bowtie2_align:
        """
        Align reads against a genome (index) with bowtie2, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bowtie2/", **config)
        output:
            get_alignment_pipes()
        log:
            expand("{log_dir}/bowtie2_align/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bowtie2_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['align']
        threads: 20
        conda:
            "../envs/bowtie2.yaml"
        shell:
            """
            bowtie2 {params.params} --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'bwa':
    config['bwaindex_types'] = ['amb', 'ann', 'bwt', 'pac', 'sa']

    rule bwa_index:
        """
        Make a genome index for bwa.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
        log:
            expand("{log_dir}/bwa_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/bwa_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            prefix="{genome_dir}/{{assembly}}/index/bwa/{{assembly}}".format(**config),
            params=config['index']
        priority: 1
        resources:
            mem_gb=5
        conda:
            "../envs/bwa.yaml"
        shell:
            "bwa index -p {params.prefix} {params.params} {input} > {log} 2>&1"


    rule bwa_mem:
        """
        Align reads against a genome (index) with bwa, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}.{bwaindex_types}", **config)
        output:
            get_alignment_pipes()
        log:
            expand("{log_dir}/bwa_mem/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bwa_mem/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            index_dir=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}", **config),
            params=config['align']
        resources:
            mem_gb=23
        threads: 20
        conda:
            "../envs/bwa.yaml"
        shell:
            """
            bwa mem {params.params} -t {threads} {params.index_dir} {input.reads} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'hisat2':
    rule hisat2_index:
        """
        Make a genome index for hisat2.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/hisat2/", **config))
        log:
            expand("{log_dir}/hisat2_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/hisat2_index/{{assembly}}.benchmark.txt", **config)[0]
        priority: 1
        threads: 4
        conda:
            "../envs/hisat2.yaml"
        params: config['index']
        shell:
            "hisat2-build {params} -p {threads} {input} {output}/{wildcards.assembly} > {log} 2>&1"


    rule hisat2_align:
        """
        Align reads against a genome (index) with hisat2, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/hisat2/", **config)
        output:
            get_alignment_pipes()
        log:
            expand("{log_dir}/hisat2_align/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/hisat2_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['align']
        threads: 20
        conda:
            "../envs/hisat2.yaml"
        shell:
            """
            hisat2 {params.params} --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'star' or config.get('quantifier', '') == 'star':
    rule star_index:
        """
        Make a genome index for STAR.
        
        Troubleshooting:
        1) sufficient disk space?
        2) increase the RAM available (--limitGenomeGenerateRAM)
        3) reduce the number of threads (snakemake -j 5)
        4) reduce accuracy (--genomeSAsparseD 2)
                
        For example, in your config.yaml, set aligner/quantifier:
        aligner:
            star:
                index: --limitGenomeGenerateRAM 60000000000 --genomeSAsparseD 1
        """
        input:
            genome = expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            sizefile= expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config),
            gtf = expand("{genome_dir}/{{assembly}}/{{assembly}}.gtf", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/{aligner}", **config))
        log:
            expand("{log_dir}/{aligner}_index/{{assembly}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/{aligner}_index/{{assembly}}.benchmark.txt", **config)[0]
        params:
            config['index']
        priority: 1
        threads: 10
        resources:
            mem_gb=37
        conda:
            "../envs/star.yaml"
        shell:
            """
            function log2 {{
                    local x=0
                    for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
                        let x=$x+1
                    done
                    echo $x
            }}
            
            # set genome dependent variables
            NBits=""
            NBases=""
            GenomeLength=$(awk -F"\t" '{{x+=$2}}END{{printf "%i", x}}' {input.sizefile})
            NumberOfReferences=$(awk 'END{{print NR}}' {input.sizefile})
            if [ $NumberOfReferences -gt 5000 ]; then
                # for large genomes, --genomeChrBinNbits should be scaled to min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])
                # ReadLength is skipped here, as it is unknown
                LpR=$(log2 $((GenomeLength / NumberOfReferences)))
                NBits="--genomeChrBinNbits $(($LpR<18 ? $LpR : 18))"
                printf "NBits: $NBits\n\n" >> {log} 2>&1
            fi
            
            if [ $GenomeLength -lt 268435456 ]; then
                # for small genomes, --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2-1)
                logG=$(( $(log2 $GenomeLength) / 2 - 1 ))
                NBases="--genomeSAindexNbases $(( $logG<14 ? $logG : 14 ))"
                printf "NBases: $NBases\n\n" >> {log} 2>&1
            fi
            
            mkdir -p {output}
            
            STAR --runMode genomeGenerate --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} \
            --genomeDir {output} --outFileNamePrefix {output}/ \
            --runThreadN {threads} $NBits $NBases {params} >> {log} 2>&1
            """


    if config.get('run_alignment', True):
        rule star_align:
            """
            Align reads against a genome (index) with STAR, and pipe the output to the required sorter(s).
            """
            input:
                reads=get_reads,
                index=expand("{genome_dir}/{{assembly}}/index/{aligner}", **config)
            output:
                dir =directory(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}", **config)),
                pipe=get_alignment_pipes()
            log:
                expand("{log_dir}/{aligner}_align/{{assembly}}-{{sample}}.log", **config)
            benchmark:
                expand("{benchmark_dir}/{aligner}_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
            params:
                input=lambda wildcards, input: f' {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                               f' {input.reads[0]} {input.reads[1]}',
                params=config['align']
            threads: 1
            resources:
                mem_gb=30
            conda:
                "../envs/star.yaml"
            shell:
                """
                mkdir -p {output.dir}
                
                STAR --genomeDir {input.index} --readFilesIn {params.input} --quantMode GeneCounts \
                --outFileNamePrefix {output.dir}/ --runThreadN {threads} {params.params} \
                --outSAMtype BAM Unsorted --outStd BAM_Unsorted | tee {output.pipe} 1> /dev/null 2>> {log}
                """


rule sambamba_sort:
    """
    Sort the result of alignment with the sambamba sorter.
    """
    input:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.sambamba-{{sorting}}.pipe", **config)
    output:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.sambamba-{{sorting}}.bam", **config)
    log:
        expand("{log_dir}/sambamba_sort/{{assembly}}-{{sample}}-sambamba_{{sorting}}.log", **config)
    group: 'alignment'
    benchmark:
        expand("{benchmark_dir}/sambamba_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        lambda wildcards: "-n" if wildcards.sorting == 'queryname' else '',
    threads: 4
    conda:
        "../envs/sambamba.yaml"
    shell:
        """
        sambamba view --nthreads {threads} -S -f bam  {input[0]} -o /dev/stdout  2> {log} |
        sambamba sort --nthreads {threads} {params}   /dev/stdin -o {output[0]}  2> {log}
        """


rule samtools_sort:
    """
    Sort the result of alignment with the samtools sorter.
    """
    input:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-{{sorting}}.pipe", **config)
    output:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-{{sorting}}.bam", **config)
    log:
        expand("{log_dir}/samtools_sort/{{assembly}}-{{sample}}-samtools_{{sorting}}.log", **config)
    group: 'alignment'
    benchmark:
        expand("{benchmark_dir}/samtools_sort/{{assembly}}-{{sample}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        order=lambda wildcards: "-n" if wildcards.sorting == 'queryname' else '',
        threads=lambda wildcards, input, output, threads: threads - 1
    threads: 4
    resources:
        mem_mb=2500
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        trap \"rm -f {output}*\" INT;
        samtools sort -@ {params.threads} {params.order} {input} -o {output}  2> {log}
        """


rule samtools_index:
    """
    Create an index of a bam file which can be used for e.g. visualization.
    """
    input:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config)
    output:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam.bai", **config)
    log:
        expand("{log_dir}/samtools_index/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/samtools_index/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        config['samtools_index']
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {params} {input} {output}
        """


rule mark_duplicates:
    """
    Mark (but keep) all duplicate reads in a bam file with picard MarkDuplicates
    """
    input:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config)
    output:
        bam=    expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config),
        metrics=expand("{qc_dir}/dedup/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.metrics.txt", **config)
    log:
        expand("{log_dir}/mark_duplicates/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/mark_duplicates/{{assembly}}-{{sample}}-{{sorter}}-{{sorting}}.benchmark.txt", **config)[0]
    params:
        config['markduplicates']
    resources:
        mem_gb=5
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates {params} INPUT={input} "
        "OUTPUT={output.bam} METRICS_FILE={output.metrics} > {log} 2>&1"
