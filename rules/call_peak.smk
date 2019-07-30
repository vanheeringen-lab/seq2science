# def aggregate_replicates_genrich(wildcards):
#     sample = samples.index[(samples['assembly']  == wildcards.assembly).values &
#                            (samples['condition'] == wildcards.condition).values &
#                            (samples['project']   == wildcards.project).values]
#     return expand(f"{{result_dir}}/{dedup_dir}/{{sample}}-{wildcards.condition}-{wildcards.project}-{wildcards.assembly}.bam",
#                   result_dir=config['result_dir'], sample=sample)

rule genrich_pileup:
    input:
        expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.bam", **config)
    output:
        bedgraphish=expand("{result_dir}/genrich/{{sample}}-{{assembly}}.bdgish", **config),
        log=expand("{result_dir}/genrich/{{sample}}-{{assembly}}.log", **config)
    log:
        expand("{log_dir}/genrich_pileup/{{sample}}-{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/genrich_pileup/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
    conda:
        "../envs/call_peak.yaml"
    params:
        config['peak_caller'].get('genrich', " ")
    threads: 15
    shell:
        "input=$(echo {input} | tr ' ' ','); "
        "Genrich -X -t $input -f {output.log} -k {output.bedgraphish} {params}"


rule call_peak_genrich:
    input:
        log=expand("{result_dir}/genrich/{{sample}}-{{assembly}}.log", **config)
    output:
        narrowpeak= expand("{result_dir}/genrich/{{sample}}-{{assembly}}_peaks.narrowPeak", **config)
    log:
        expand("{log_dir}/call_peak_genrich/{{sample}}-{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/call_peak_genrich/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
    conda:
        "../envs/call_peak.yaml"
    params:
        config['peak_caller']['genrich']
    threads: 1
    shell:
        "Genrich -P -f {input.log} -o {output.narrowpeak} {params}"


config['macs2_types'] = ['control_lambda.bdg', 'summits.bed', 'peaks.narrowPeak',
                         'peaks.xls', 'treat_pileup.bdg']
def get_fastqc(wildcards):
    if config['layout'][wildcards.sample] == "SINGLE":
        return expand("{result_dir}/{trimmed_dir}/SE/{{sample}}_trimmed_fastqc.zip", **config)
    return sorted(expand("{result_dir}/{trimmed_dir}/PE/{{sample}}_{fqext1}_trimmed_fastqc.zip", **config))

rule call_peak_macs2:
    #
    # Calculates genome size based on unique kmers of average length
    #
    input:
        bam=   expand("{result_dir}/{dedup_dir}/{{sample}}-{{assembly}}.bam", **config),
        fastqc=get_fastqc
    output:
        expand("{result_dir}/macs2/{{sample}}-{{assembly}}_{macs2_types}", **config)
    log:
        expand("{log_dir}/call_peak_macs2/{{sample}}-{{assembly}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/call_peak_macs2/{{sample}}-{{assembly}}.benchmark.txt", **config)[0]
    params:
        name=lambda wildcards, input: f"{wildcards.sample}" if config['layout'][wildcards.sample] == 'SINGLE' else \
                                      f"{wildcards.sample}_{config['fqext1']}",
        genome=f"{config['genome_dir']}/{{assembly}}/{{assembly}}.fa",
        macs_params=config['peak_caller']['macs2']
    conda:
        "../envs/call_peak_macs2.yaml"
    shell:
        f"""
        # extract the kmer size, and get the effective genome size from it
        kmer_size=$(unzip -p {{input.fastqc}} {{params.name}}_trimmed_fastqc/fastqc_data.txt  | grep -P -o '(?<=Sequence length\\t).*' | grep -P -o '\d+$');
        GENSIZE=$(unique-kmers.py {{params.genome}} -k $kmer_size --quiet 2>&1 | grep -P -o '(?<=\.fa: ).*');
        echo "kmer size: $kmer_size, and effective genome size: $GENSIZE" >> {{log}}

        # call peaks
        macs2 callpeak -t {{input.bam}} --outdir {config['result_dir']}/macs2/ -n {{wildcards.sample}}-{{wildcards.assembly}} \
        {{params.macs_params}} -g $GENSIZE > {{log}} 2>&1
        """
