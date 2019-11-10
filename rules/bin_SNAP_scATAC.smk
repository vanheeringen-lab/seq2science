rule create_bins_SNAP_object:
    '''Add a Binned genome matrix with 5kb bins to the SNAPobject, after which it is renamed and moved to the Snapfiles folder for downstream analysis in R using Snaptools'''
    input:
        config['result_dir'] + '{plate}/merged_ordered.snap'
    output:
        config['result_dir'] + 'Snapfiles/{plate}_binned_ordered.snap'
    conda:
        "../envs/Snaptools.yml"
    params:
        bin_opt=config['bin_opt']
    shell:
        ''' 
        snaptools snap-add-bmat --snap-file={input} {params.bin_opt}
        echo 'bmat added,moving file'
        mv {input} {output}	
        '''