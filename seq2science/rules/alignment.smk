def get_reads(wildcards):
    """
    Function that returns the reads for any aligner.
    """
    if 'replicate' in samples and config.get('technical_replicates') == 'merge':
        if config['layout'].get(wildcards.sample, False) == "SINGLE":
            return expand("{trimmed_dir}/merged/{{sample}}_trimmed.{fqsuffix}.gz", **config)
        return sorted(expand("{trimmed_dir}/merged/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))
    else:
        if config['layout'].get(wildcards.sample, False) == "SINGLE":
            return expand("{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config)
        return sorted(expand("{trimmed_dir}/{{sample}}_{fqext}_trimmed.{fqsuffix}.gz", **config))


if config['aligner'] == 'bowtie2':
    rule bowtie2_index:
        """
        Make a genome index for bowtie2. This index is required for faster alignment.
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
            pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])
        log:
            expand("{log_dir}/bowtie2_align/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bowtie2_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['align']
        threads: 9
        conda:
            "../envs/bowtie2.yaml"
        shell:
            """
            bowtie2 {params.params} --threads {threads} -x {input.index}{wildcards.assembly} {params.input} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'bwa-mem':
    rule bwa_index:
        """
        Make a genome index for bwa (mem). This index is required for faster alignment.
        """
        input:
            expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config)
        output:
            directory(expand("{genome_dir}/{{assembly}}/index/bwa/", **config)),
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
        Align reads against a genome (index) with bwa-mem, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/bwa/", **config)
        output:
            pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])
        log:
            expand("{log_dir}/bwa_mem/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/bwa_mem/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            index_dir=expand("{genome_dir}/{{assembly}}/index/bwa/{{assembly}}", **config),
            params=config['align']
        resources:
            mem_gb=13
        threads: 10
        conda:
            "../envs/bwa.yaml"
        shell:
            """
            bwa mem {params.params} -t {threads} {params.index_dir} {input.reads} 2> {log} | tee {output} 1> /dev/null 2>> {log}
            """


elif config['aligner'] == 'hisat2':
    rule hisat2_index:
        """
        Make a genome index for hisat2. This index is required for faster alignment.
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
            pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])
        log:
            expand("{log_dir}/hisat2_align/{{assembly}}-{{sample}}.log", **config)
        group: 'alignment'
        benchmark:
            expand("{benchmark_dir}/hisat2_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f'-U {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f'-1 {input.reads[0]} -2 {input.reads[1]}',
            params=config['align']
        threads: 9
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
            gtf = expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config)
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


    rule star_align:
        """
        Align reads against a genome (index) with STAR, and pipe the output to the required sorter(s).
        """
        input:
            reads=get_reads,
            index=expand("{genome_dir}/{{assembly}}/index/{aligner}", **config)
        output:
            dir =directory(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}", **config)),
            pipe=pipe(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)[0])
        log:
            directory(expand("{log_dir}/{aligner}_align/{{assembly}}-{{sample}}", **config))
        benchmark:
            expand("{benchmark_dir}/{aligner}_align/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
        params:
            input=lambda wildcards, input: f' {input.reads}' if config['layout'][wildcards.sample] == 'SINGLE' else \
                                           f' {input.reads[0]} {input.reads[1]}',
            params=config['align']
        threads: 8
        resources:
            mem_gb=30
        conda:
            "../envs/star.yaml"
        shell:
            """
            trap "find {log} -type f ! -name Log* -exec rm {{}} \;" EXIT
            mkdir -p {log}
            mkdir -p {output.dir}                

            STAR --genomeDir {input.index} --readFilesIn {params.input} --readFilesCommand gunzip -c \
            --quantMode GeneCounts --outSAMtype BAM Unsorted --outStd BAM_Unsorted \
            --outFileNamePrefix {log}/ --outTmpDir {output.dir}/STARtmp \
            --runThreadN {threads} {params.params} > {output.pipe} 2> {log}/Log.stderr.out

            # move all non-log files to output directory (this way the log files are kept on error)
            find {log} -type f ! -name Log* -exec mv {{}} {output.dir} \;
            """


rule samtools_presort:
    """
    (Pre)sort the result of alignment with the samtools sorter.
    """
    input:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.pipe", **config)
    output:
        temp(expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate-unsieved.bam", **config))
    log:
        expand("{log_dir}/samtools_presort/{{assembly}}-{{sample}}.log", **config)
    group: 'alignment'
    benchmark:
        expand("{benchmark_dir}/samtools_presort/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    params:
        threads=lambda wildcards, input, output, threads: max([1, threads - 1]),
        sort_order="-n" if not sieve_bam(config) and config.get("bam_sort_order") == "queryname" else "",
        out_dir=f"{config['result_dir']}/{config['aligner']}"
    threads: 2
    resources:
        mem_mb=2500
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        trap "rm -f {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp*bam" INT;
        samtools sort {params.sort_order} -@ {params.threads} {input} -o {output} -T {params.out_dir}/{wildcards.assembly}-{wildcards.sample}.tmp 2> {log}
        """
