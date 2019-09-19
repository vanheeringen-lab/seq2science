def get_strandedness(wildcards):
    sample = f"{wildcards.sample}"
    strandedness = samples["strandedness"].loc[sample]
    return strandedness

rule bam_stranded_bigwig:
    """
    Convert a bam file into two bigwig files, one for each strand    
    """
    input:
        bam=expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config),
        bai=expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bai", **config)
    output:
        forward=expand("{result_dir}/bigwigs/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.fwd.bw", **config),
        reverse=expand("{result_dir}/bigwigs/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.rev.bw", **config),
    params:
        flags=config['bam_bigwig']['deeptools'],
        strandedness=get_strandedness
    wildcard_constraints:
        sample=any_sample(),
        sorting=config['bam_sort_order']
    log:
        expand("{log_dir}/bam_bigwig/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bam_bigwig/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.log", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 20
    resources:
        deeptools_limit=1
    shell:
        """       
        direction1=forward
        direction2=reverse
        if [ {params.strandedness} == 'reverse' ]; then
            direction1=reverse
            direction2=forward
        fi
                    
        bamCoverage --bam {input.bam} --outFileName {output.forward} --filterRNAstrand $direction1 --numberOfProcessors {threads} {params.flags} --verbose >> {log} 2>&1 &&        
        bamCoverage --bam {input.bam} --outFileName {output.reverse} --filterRNAstrand $direction2 --numberOfProcessors {threads} {params.flags} --verbose >> {log} 2>&1
        """

rule bam_bigwig:
    """
    Convert a bam file into a bigwig file
    """
    input:
        expand("{dedup_dir}/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bam", **config)
    output:
        expand("{result_dir}/bigwigs/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.bw", **config)
    params:
        config['bam_bigwig']['deeptools']
    wildcard_constraints:
        sample=any_sample(),
        sorting=config['bam_sort_order']
    log:
        expand("{log_dir}/bam_bigwig/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/bam_bigwig/{{assembly}}-{{sample}}.{{sorter}}-{{sorting}}.log", **config)[0]
    conda:
        "../envs/deeptools.yaml"
    threads: 20
    resources:
        deeptools_limit=1
    shell:
        """
        bamCoverage --bam {input} --outFileName {output} --numberOfProcessors {threads} {params} --verbose >> {log} 2>&1
        """
