rule trim_galore_PE:
    '''
    Automated adapter detection, adapter trimming, and quality trimming through trim galore (paired-end).
    '''	
    input:
        r1= config['raw_data_loc'] + '{plate}/{cell}_R1.fastq.gz',
        r2= config['raw_data_loc'] + '{plate}/{cell}_R2.fastq.gz'
    output:
        r1= config['result_dir'] + '{plate}/trim_fq/{cell}_R1_val_1.fq.gz',
        r2= config['result_dir'] + '{plate}/trim_fq/{cell}_R2_val_2.fq.gz'
    conda:
        "../envs/trimgalore.yaml"
    params:
        trim_galore_settings = config['trim_galore']
    shell:
        '''trim_galore {params.trim_galore_settings}\
         -o $(dirname {output.r1}) \
        {input.r1} {input.r2}'''