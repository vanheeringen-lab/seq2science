rule sort_plates_bam:
    '''Sort the plate bamfiles'''
    input:
        output_loc + '{plate}/f2q30_merged.bam'
    output:
        output_loc + '{plate}/f2q30_merged_sorted.bam'
    conda:
        "envs/python3_2.yml"
    shell:
        '''samtools sort -n {input}>{output}
