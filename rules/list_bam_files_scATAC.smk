rule list_bam:
    '''create a txt file containing a list of all BAM files per plate to merge the cingle cell BAM into a whole plate BAM'''
    params:
        samples_p = samples['plate'],
	samples_c = samples['cell'],
        plate_vallues  = plate_vallues
    input:
        expand(config['result_dir'] + '{plate}/cell_ID_BAMs/{cell}.bam', zip,
               plate= {params.samples_p},
               cell= {params.samples_c})
    output:
        expand(config['result_dir'] + '{plate}/bam_file_list.txt', plate= {params.plate_vallues})
    conda:
        "../envs/python3.yml"
    shell:
        ''' ../scripts/list_bam.sh {output_loc}
        '''
