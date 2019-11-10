rule sort_plates_bam:
    '''Sort the plate bamfiles'''
    input:
        config['result_dir'] + '{plate}/f2q30_merged.bam'
    output:
        config['result_dir'] + '{plate}/f2q30_merged_sorted.bam'
    conda:
        "../envs/python3_2.yml"
    shell:
        '''samtools sort -n {input}>{output}'''
