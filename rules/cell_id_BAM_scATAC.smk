rule cell_id_BAM:
    '''add the cell ID to the Qname of the BAMfile of each cell'''
    input:
        config['result_dir'] + '{plate}/hisat2_hg38_mapped/{cell}_f2q30.bam',
    output:
        config['result_dir'] + '{plate}/cell_ID_BAMs/{cell}.bam'
    conda:
        "../envs/python3.yml"
    script:
        '../scripts/add_cell_ID.py'