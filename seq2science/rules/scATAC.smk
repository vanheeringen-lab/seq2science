"""
all rules/logic specific to the single cell ATAC workflow should be here.
"""


rule create_SNAP_object:
    """
    Create a snapobject for each BAM file. 

    These snapobjects can be merged later using snaptools in R.
    """
    input:
        bams=expand("{final_bam_dir}/{{assembly}}-{{sample}}.sambamba-queryname.bam", **config),
        genome_size=rules.get_genome_support_files.output.sizes,
    output:
        expand("{result_dir}/snap/{{assembly}}-{{sample}}.snap", **config),
    log:
        expand("{log_dir}/create_SNAP_object/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/create_SNAP_object/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    threads: 4
    conda:
        "../envs/snaptools.yaml"
    params:
        params=config["snaptools_opt"],
        chrm=f"--keep-chrm={'TRUE' if not config['remove_mito'] else 'FALSE'}",
        mapq=f"--min-mapq={config['min_mapping_quality']}",
    shell:
        """
        snaptools snap-pre --input-file={input.bams} --output-snap={output} --genome-name={wildcards.assembly} \
        --genome-size={input.genome_size} {params.params} {params.chrm} {params.mapq} > {log} 2>&1
        """


rule create_bins_SNAP_object:
    """
    Add a Binned genome matrix to the SNAPobject, after which it is renamed and moved
    to the Snapfiles folder for downstream analysis in R using Snaptools
    """
    input:
        rules.create_SNAP_object.output,
    output:
        expand("{result_dir}/snap/{{assembly}}-{{sample}}.binned.snap", **config),
    log:
        expand("{log_dir}/create_bins_SNAP_object/{{assembly}}-{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/create_SNAP_object/{{assembly}}-{{sample}}.benchmark.txt", **config)[0]
    conda:
        "../envs/snaptools.yaml"
    message: EXPLAIN["create_bins_SNAP_object"]
    params:
        config["bin_opt"],
    shell:
        """ 
        snaptools snap-add-bmat --snap-file={input} {params} > {log} 2>&1
        echo "bmat added, moving file" >> {log}
        mv {input} {output} >> {log} 2>&1
        """
