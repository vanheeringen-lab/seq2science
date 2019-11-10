rule merge_plates:
    '''Merge the bam files of individual cells (containing the cell ID in the qName of the BAM file) into one bam file per plate'''
    input:
        config['result_dir'] + '{plate}/bam_file_list.txt'
    output:
        config['result_dir'] + '{plate}/f2q30_merged.bam'
    conda:
        "../envs/python3.yml"
    shell:
        ''' samtools merge -b {input} {output}
        ''' 