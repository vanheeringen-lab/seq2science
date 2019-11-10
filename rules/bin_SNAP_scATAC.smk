rule create_bins_SNAP_object:
    '''Add a Binned genome matrix with 5kb bins to the SNAPobject, after which it is renamed and moved to the Snapfiles folder for downstream analysis in R using Snaptools'''
    input:
        output_loc + '{plate}/merged_ordered.snap'
    output:
        output_loc + 'Snapfiles/{plate}_binned_ordered.snap'
    conda:
        "envs/Snaptools.yml"
    shell:
        ''' 
        snaptools snap-add-bmat --snap-file={input} --bin-size-list 5000 --verbose=True
        echo 'bmat added,moving file'
        mv {input} {output}	
        '''