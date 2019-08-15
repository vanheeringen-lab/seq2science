def get_genrich_replicates(wildcards):
    sample_condition, assembly = wildcards.fname.split('-')
    if not 'condition' in samples.columns \
    or (sample_condition in samples.index and not sample_condition in samples['condition'].values):
        return expand(f"{{result_dir}}/{{dedup_dir}}/{wildcards.fname}.bam", **config)
    else:
        return expand([f"{{result_dir}}/{{dedup_dir}}/{replicate}-{assembly}.bam"
        for replicate in samples[(samples['assembly'] == assembly) & (samples['condition'] == sample_condition)].index], **config)


rule genrich_pileup:
    input:
        get_genrich_replicates
    output:
        bedgraphish=expand("{result_dir}/genrich/{{fname}}.bdgish", **config),
        log=expand("{result_dir}/genrich/{{fname}}.log", **config)
    log:
        expand("{log_dir}/genrich_pileup/{{fname}}_pileup.log", **config)
    benchmark:
        expand("{benchmark_dir}/genrich_pileup/{{fname}}.benchmark.txt", **config)[0]
    conda:
        "../envs/genrich.yaml"
    params:
        config['peak_caller'].get('genrich', " ")  # TODO: move this to config.schema.yaml
    threads: 15  # TODO: genrich uses lots of ram. Get the number from benchmark, instead of doing it through threads
    shell:
        """
        input=$(echo {input} | tr ' ' ',')
        Genrich -X -t $input -f {output.log} -k {output.bedgraphish} {params} -v > {log} 2>&1
        """


rule call_peak_genrich:
    input:
        log=expand("{result_dir}/genrich/{{fname}}.log", **config)
    output:
        narrowpeak= expand("{result_dir}/genrich/{{fname}}_peaks.narrowPeak", **config)
    log:
        expand("{log_dir}/call_peak_genrich/{{fname}}_peak.log", **config)
    benchmark:
        expand("{benchmark_dir}/call_peak_genrich/{{fname}}.benchmark.txt", **config)[0]
    conda:
        "../envs/genrich.yaml"
    params:
        config['peak_caller'].get('genrich', "")
    threads: 1
    shell:
        "Genrich -P -f {input.log} -o {output.narrowpeak} {params} -v > {log} 2>&1"


config['macs2_types'] = ['control_lambda.bdg', 'summits.bed', 'peaks.narrowPeak',
                         'peaks.xls', 'treat_pileup.bdg']
def get_fastqc(wildcards):
    if config['layout'][wildcards.sample] == "SINGLE":
        return expand("{result_dir}/{trimmed_dir}/{{sample}}_trimmed_fastqc.zip", **config)
    return sorted(expand("{result_dir}/{trimmed_dir}/{{sample}}_{fqext1}_trimmed_fastqc.zip", **config))

rule macs2_callpeak:
    #
    # Calculates genome size based on unique kmers of average length
    #
    input:
        bam=expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.bam", **config),
        fastqc=get_fastqc
    output:
        expand("{result_dir}/macs2/{{sample}}-{{assembly}}_{macs2_types}", **config)
    log:
        expand("{log_dir}/macs2_callpeak/{{sample}}-{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/macs2_callpeak/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
    params:
        name=lambda wildcards, input: f"{wildcards.sample}" if config['layout'][wildcards.sample] == 'SINGLE' else \
                                      f"{wildcards.sample}_{config['fqext1']}",
        genome=f"{config['genome_dir']}/{{assembly}}/{{assembly}}.fa",
        macs_params=config['peak_caller'].get('macs2', "")  # TODO: move to config.schema.yaml
    conda:
        "../envs/macs2.yaml"
    shell:
        f"""
        # extract the kmer size, and get the effective genome size from it
        kmer_size=$(unzip -p {{input.fastqc}} {{params.name}}_trimmed_fastqc/fastqc_data.txt  | grep -P -o '(?<=Sequence length\\t).*' | grep -P -o '\d+$');
        GENSIZE=$(unique-kmers.py {{params.genome}} -k $kmer_size --quiet 2>&1 | grep -P -o '(?<=\.fa: ).*');
        echo "kmer size: $kmer_size, and effective genome size: $GENSIZE" >> {{log}}

        # call peaks
        macs2 callpeak --bdg -t {{input.bam}} --outdir {config['result_dir']}/macs2/ -n {{wildcards.sample}}-{{wildcards.assembly}} \
        {{params.macs_params}} -g $GENSIZE > {{log}} 2>&1
        """
