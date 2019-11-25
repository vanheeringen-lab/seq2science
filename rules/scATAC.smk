def all_cells_for_plate(wildcards):
    import glob
    return glob.glob(expand(f"{{fastq_dir}}/{wildcards.sample}/*{{fqsuffix}}.gz", **config)[0])


rule cell_ID2_fastq_ID:
    '''
    takes the final part of the file name, between the - and _ characters and uses this as a cell identifier.
    '''
    input:
        all_cells_for_plate
    output:
        R1=temp(expand("{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config)),
        R2=temp(expand("{fastq_dir}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config))
    run:
        cells = set(cell[:-(3+len(config['fqsuffix'])+1+len(config['fqext1']))] for cell in input)
        assert len(cells) * 2 == len(input)
        for cell in list(cells):
            cell_f = f"{cell}{config['fqext1']}.{config['fqsuffix']}.gz"
            cell_r = f"{cell}{config['fqext2']}.{config['fqsuffix']}.gz"
            # once for forward
            shell(f"""cell_id='@'$(echo {cell_f} | awk -F"_" '{{{{ print $2 }}}}'| awk -F"-" '{{{{ print $NF }}}}')':'; """
                  f"""zcat {cell_f}|awk -v var=$cell_id '{{{{gsub(/@/, var )}}}}1' | gzip >> {output.R1}""")
            # once for forward
            shell(f"""cell_id='@'$(echo {cell_r} | awk -F"_" '{{{{ print $2 }}}}'| awk -F"-" '{{{{ print $NF }}}}')':'; """
                  f"""zcat {cell_r}|awk -v var=$cell_id '{{{{gsub(/@/, var )}}}}1' | gzip >> {output.R2}""")


# rule cell_id_BAM:
#     '''
#     Add the cell ID to the Qname of the BAMfile of each cell.
#     '''
#     input:
#         expand("{dedup_dir}/{{assembly}}-{{sample}}.sambamba-queryname.bam", **config)
#     output:
#         expand("{result_dir}/cell_id/{{assembly}}-{{sample}}.sambamba-queryname.bam", **config)
#     run:
#         import pysam
#         import os
#
#         path = '/'.join(output[0].split('/')[:-1])
#         if not os.path.exists(path):
#             os.makedirs(path)
#
#         samfile = pysam.AlignmentFile(input[0], "rb")
#         sam_out = pysam.AlignmentFile(output[0], "wb", template=samfile)
#
#         for read in samfile:
#             read.qname = f"{wildcards.sample}-{read.qname}"
#             sam_out.write(read)
#
#         sam_out.close()
#         samfile.close()
#
#
# def get_all_cells_per_plate(wildcards):
#     """
#     Function that returns all bams that belong to a plate
#     """
#     cells = [cell for cell in samples.index if wildcards.plate in cell]
#     assert len(set(samples['assembly'])) == 1, "this logic has not been implemented (yet) for " \
#                                                "multiple assemblies at the same time."
#
#     bams = [expand(f"{{result_dir}}/cell_id/{wildcards.assembly}-{cell}.sambamba-queryname.bam", **config)[0] for cell in cells]
#     return bams
#
#
# rule merge_plates:
#     '''
#     Merge the bam files of individual cells (containing the cell ID in the qName of the BAM file)
#     into one bam file per plate
#     '''
#     input:
#         get_all_cells_per_plate
#     output:
#         expand("{dedup_dir}/{{assembly}}-{{plate,plate\d}}.merged.sambamba-queryname.bam", **config)
#     log:
#         expand("{log_dir}/merge_plates/{{assembly}}-{{plate}}.log", **config)
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         '''
#         samtools merge {output} {input} > {log} 2>&1
#         '''
#
#
# def get_plate_bam(wildcards):
#     assert len(set(samples['assembly'])) == 1, "this logic has not been implemented (yet) for " \
#                                                "multiple assemblies at the same time."
#     return expand(f"{{dedup_dir}}/{samples['assembly'][0]}-{wildcards.plate}.merged.sambamba-queryname.bam", **config)


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
