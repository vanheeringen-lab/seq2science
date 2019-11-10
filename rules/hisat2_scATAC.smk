rule hisat2:
    '''Use hisat2 to map the reads files to the genome hg38'''
    input:
        r1= config['result_dir'] + '{plate}/trim_fq/{cell}_R1_val_1.fq.gz',
        r2= config['result_dir'] + '{plate}/trim_fq/{cell}_R2_val_2.fq.gz'
    output:
        bam= config['result_dir'] + '{plate}/hisat2_hg38_mapped/{cell}_f2q30.bam',
        stats= config['result_dir'] + '{plate}/hisat2_hg38_log/{cell}_aln_sum.txt'
    threads: 4
    params:
        genome= config['genome_dir'],
        hisat_opt = config['hisat2'],
        samtools_options = config['samtools']
    conda:
        "../envs/hisat2.yml"
    shell:
        ''' 
        hisat2 {params.hisat_opt} -p {threads} -x {params.genome} -1 {input.r1} -2 {input.r2} --summary-file {output.stats} | \
            samtools view {params.samtools_options} | \
            samtools sort - -T {wildcards.cell}_tmp -o {output.bam}
        '''