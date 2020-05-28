rule create_SNAP_object:
    """
    Create a snapobject for each BAM file. 
    
    These snapobjects can be merged later using snaptools in R.
    """
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
        params=config["snaptools_opt"],
        assembly=samples["assembly"][0]
    shell:
        """
        snaptools snap-pre --input-file={input.bams} --output-snap={output} --genome-name={params.assembly} \
        --genome-size={input.genome_size} {params.params} > {log} 2>&1
        """

rule create_bins_SNAP_object:
    """
    Add a Binned genome matrix with 5kb bins to the SNAPobject, after which it is renamed and moved
    to the Snapfiles folder for downstream analysis in R using Snaptools
    """
    input:
        expand("{result_dir}/snap/{{sample}}-{{assembly}}.snap", **config)
    output:
        expand("{result_dir}/snap/{{sample}}-{{assembly}}.binned.snap", **config)
    log:
        expand("{log_dir}/create_bins_SNAP_object/{{sample}}-{{assembly}}.log", **config)
    conda:
        "../envs/snaptools.yaml"
    params:
        config["bin_opt"]
    shell:
        """ 
        snaptools snap-add-bmat --snap-file={input} {params} > {log} 2>&1
        echo "bmat added, moving file" > {log} 2>&1
        mv {input} {output}	> {log} 2>&1
        """
