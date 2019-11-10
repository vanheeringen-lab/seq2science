from glob import iglob
import pandas as pd

#I make a list of all folders (./plate1, ./plate2 etc.) containing a fastq folder containing 'cell_identifier.fastq.gz' 
r1 = iglob(str(config['raw_data_loc'] + '*/*_R1.fastq.gz'))
samples = pd.DataFrame()
samples['r1'] = [i for i in r1]

#I save the variable of the folder name as a list of 'plate' variables, while the cell_identifier is saved as a list of 'cell' variables. These variables will be used in snakemake to make sure for each cell the pipeline runs correctly
tmp = samples["r1"].str.split("/", expand = True) 
samples[['plate', 'cell']] = tmp[tmp.columns[-2:]]
samples['cell'].replace(regex=True,inplace=True,to_replace=r'_R1.fastq.gz',value=r'')
plate_vallues  = list(set(samples['plate'].values))

rule list_bam:
    '''create a txt file containing a list of all BAM files per plate to merge the cingle cell BAM into a whole plate BAM'''
    params:
        outdir = config['result_dir']
    conda:
        "../envs/python3.yml"
    input:
        expand(config['result_dir'] + '{plate}/cell_ID_BAMs/{cell}.bam', zip,
               plate=samples['plate'],
               cell=samples['cell'])
    output:
        expand(config['result_dir'] + '{plate}/bam_file_list.txt', plate= plate_vallues)
    shell:
        ''' ../../scripts/list_bam.sh {params.outdir}
        '''