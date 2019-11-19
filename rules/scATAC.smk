rule cell_id_BAM:
    '''
    Add the cell ID to the Qname of the BAMfile of each cell.
    '''
    input:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.bam", **config)
    output:
        expand("{result_dir}/{aligner}/{{assembly}}-{{sample}}.samtools-coordinate.cell_id.bam", **config)
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
    assert len(set(samples['assembly'])) == 1, "this logic has not been implemented (yet) for " \
                                               "multiple assemblies at the same time."

    bams = [expand(f"{{result_dir}}/{{aligner}}/{wildcards.assembly}-{cell}.samtools-coordinate.cell_id.bam", **config)[0] for cell in cells]
    return bams


rule merge_plates:
    '''
    Merge the bam files of individual cells (containing the cell ID in the qName of the BAM file)
    into one bam file per plate
    '''
    input:
        get_all_cells_per_plate
    output:
        expand("{{result_dir}}/{{aligner}}/{{assembly}}-{{plate,plate\d}}.samtools-coordinate.pipe", **config)
    log:
        expand("{log_dir}/merge_plates/{{assembly}}-{{plate}}.log", **config)
    conda:
        "envs/samtools.yml"
    shell:
        ''' 
        samtools merge -b {input} {output} > {log} 2>&1
        '''


def get_plate_bam(wildcards):
    assert len(set(samples['assembly'])) == 1, "this logic has not been implemented (yet) for " \
                                               "multiple assemblies at the same time."
    return expand(f"{{dedup_dir}}/{samples['assembly'][0]}-{wildcards.plate}.samtools-coordinate.bam", **config)


rule create_SNAP_object:
    '''
    Create a snapobject for each BAM file. 
    
    These snapobjects can be merged later using snaptools in R.
    '''
    input:
        get_plate_bam
    output:
        expand("{result_dir}/snap/{{plate}}.snap", **config)
    log:
        expand("{log_dir}/create_SNAP_object/{{plate}}.log", **config)
    threads: 4
    conda:
        "../envs/Snaptools.yml"
    params:
        config['snaptools_opt']
    shell:
        '''
        snaptools snap-pre --input-file={input} --output-snap={output} {params} > {log} 2>&1
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
    log:
        expand("{log_dir}/create_bins_SNAP_object/{{plate}}.log", **config)
    conda:
        "../envs/Snaptools.yml"
    params:
        config['bin_opt']
    shell:
        ''' 
        snaptools snap-add-bmat --snap-file={input} {params} > {log} 2>&1
        echo 'bmat added, moving file' > {log} 2>&1
        mv {input} {output}	> {log} 2>&1
        '''
