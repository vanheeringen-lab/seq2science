rule create_SNAP_object:
    '''Create a snapobject for each BAM file, these snapobjects can be merged later using snaptools in R'''
    input:
        config['result_dir'] + '{plate}/f2q30_merged_sorted.bam'
    output:
        config['result_dir'] + '{plate}/merged_ordered.snap'
    threads: 4
    conda:
        "../envs/Snaptools.yml"
    params:
        snaptools_opt= config['snaptools_opt']
    shell:
        '''
        snaptools snap-pre  \
        --input-file={input}  \
        --output-snap={output}  \
        {params.snaptools_opt}
                '''