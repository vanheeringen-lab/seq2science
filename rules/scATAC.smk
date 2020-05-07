rule annotation_scATAC:
    '''Generates a file mapping all plates and wells to primer numbers. Needed for downstream cell type linking in the SNAPobject'''
    input:
        config['annotation_file']
    output:
        config['result_dir'] + '/plates_overview.csv'
    conda:
        "../envs/annotation_scATAC.yml"
    params:
        config['result_dir'] + '/scATAC_plate_primer_well.tsv'
    script:
        '../scripts/annotate_scATAC.py'


# #TODO: @Jos, is this still necessary after seqadmin demultiplexing fixes?
# rule rmv_identical_readname_artefacts_fastq_ID:
#     '''
#     Checks the readname of all fastqfiles, and removes reads with identical read names (artifact, should not be present in theory).
#     '''
#     input:
#         R1=expand("{fastq_dir}/pre-{{sample}}_{fqext1}.{fqsuffix}.gz", **config),
#         R2=expand("{fastq_dir}/pre-{{sample}}_{fqext2}.{fqsuffix}.gz", **config)
#     output:
#         R1=temp(expand("{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config)),
#         R2=temp(expand("{fastq_dir}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config))
#     conda:
#         "../envs/rmv_dup_readnames_fastq.yml"
#     log:
#         R1=expand("{log_dir}/remove_identical_fastq_Names/{{sample}}_{fqext1}.log", **config),
#         R2=expand("{log_dir}/remove_identical_fastq_Names/{{sample}}_{fqext2}.log", **config)
#     shell:
#         '''
#         seqkit rmdup {input.R1} -n -o {output.R1} > {log.R1}
#         seqkit rmdup {input.R2} -n -o {output.R2} > {log.R2}
#         '''

rule create_SNAP_object:
    '''
    Create a snapobject for each BAM file. 
    
    These snapobjects can be merged later using snaptools in R.
    '''
    input:
        bams=expand("{dedup_dir}/{{assembly}}-{{sample}}.sambamba-queryname.bam", **config),
        genome_size=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        expand("{result_dir}/snap/{{sample}}-{{assembly}}.snap", **config)
    log:
        expand("{log_dir}/create_SNAP_object/{{sample}}-{{assembly}}.log", **config)
    threads: 4
    conda:
        "../envs/snaptools.yaml"
    params:
        params=config['snaptools_opt'],
        assembly=samples['assembly'][0]
    shell:
        '''
        snaptools snap-pre --input-file={input.bams} --output-snap={output} --genome-name={params.assembly} \
        --genome-size={input.genome_size} {params.params} > {log} 2>&1
        '''

rule create_bins_SNAP_object:
    '''
    Add a Binned genome matrix with 5kb bins to the SNAPobject, after which it is renamed and moved
    to the Snapfiles folder for downstream analysis in R using Snaptools
    '''
    input:
        expand("{result_dir}/snap/{{sample}}-{{assembly}}.snap", **config)
    output:
        expand("{result_dir}/snap/{{sample}}-{{assembly}}.binned.snap", **config)
    log:
        expand("{log_dir}/create_bins_SNAP_object/{{sample}}-{{assembly}}.log", **config)
    conda:
        "../envs/snaptools.yaml"
    params:
        config['bin_opt']
    shell:
        ''' 
        snaptools snap-add-bmat --snap-file={input} {params} > {log} 2>&1
        echo 'bmat added, moving file' > {log} 2>&1
        mv {input} {output}	> {log} 2>&1
        '''
