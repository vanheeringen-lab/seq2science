rule cell_id_BAM:
    '''
    Add the cell ID to the Qname of the BAMfile of each cell.
    '''
    input:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    output:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.samtools-coordinate.cell_id.bam", **config)
    run:
        import pysam
        import os

        path = output.split('/')[:-1]
        if not os.path.exists(path):
            os.makedirs(path)

        samfile = pysam.AlignmentFile(input, "rb")
        sam_out = pysam.AlignmentFile(output, "wb", template=samfile)

        for read in samfile:
            read.qname = f"{wildcards.sample}:{read.qname}"
            sam_out.write(read)

        sam_out.close()
        samfile.close()


def get_all_cells_per_plate(wildcards):
    """
    Function that returns all bams that belong to a plate
    """
    cells = [cell for cell in samples.index if wildcards.plate in cell]
    assembly = samples['assembly'][0]
    assert len(set(samples['assembly'])) == 1, "this logic has not been implemented (yet) for " \
                                            "multiple assemblies."

    bams = [expand(f"{{dedup_dir}}/{assembly}-{cell}.samtools-coordinate.cell_id.bam", **config)[0] for cell in cells]
    return bams


rule merge_plates:
    '''
    Merge the bam files of individual cells (containing the cell ID in the qName of the BAM file)
    into one bam file per plate
    '''
    input:
        get_all_cells_per_plate
    output:
        expand("{dedup_dir}/{{plate,plate\d}}.bam", **config)
    conda:
        "envs/samtools.yml"
    shell:
        ''' 
        samtools merge -b {input} {output}
        '''

# TODO: implement sorting, best to reuse existing rule.
# TODO: not sure how this would be done most easily (we could rename {{plate}}.bam -> {{plate}}.pipe
# TODO: not sure how misleading that is :)

rule create_SNAP_object:
    '''
    Create a snapobject for each BAM file. 
    
    These snapobjects can be merged later using snaptools in R.
    '''
    input:
        expand("{dedup_dir}/{{plate}}.bam", **config)
    output:
        expand("{result_dir}/snap/{{plate}}.snap", **config)
    threads: 4
    conda:
        "../envs/Snaptools.yml"
    params:
        config['snaptools_opt']
    shell:
        '''
        snaptools snap-pre --input-file={input} --output-snap={output} {params}
        '''


rule create_bins_SNAP_object:
    '''
    Add a Binned genome matrix with 5kb bins to the SNAPobject, after which it is renamed and moved
    to the Snapfiles folder for downstream analysis in R using Snaptools
    '''
    input:
        expand("{result_dir}/snap/{{plate}}.snap", **config)
    output:
        expand("{result_dir}/snap/{{plate}}.binned.snap", **config)
    conda:
        "../envs/Snaptools.yml"
    params:
        config['bin_opt']
    shell:
        ''' 
        snaptools snap-add-bmat --snap-file={input} {params}
        echo 'bmat added, moving file'
        mv {input} {output}	
        '''
