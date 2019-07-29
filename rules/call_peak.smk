# def aggregate_replicates_genrich(wildcards):
#     sample = samples.index[(samples['assembly']  == wildcards.assembly).values &
#                            (samples['condition'] == wildcards.condition).values &
#                            (samples['project']   == wildcards.project).values]
#     return expand(f"{{result_dir}}/dedup/{{sample}}-{wildcards.condition}-{wildcards.project}-{wildcards.assembly}.bam",
#                   result_dir=config['result_dir'], sample=sample)

rule genrich_pileup:
    # TODO log
    input:
        expand("{result_dir}/dedup/{{sample}}-{{assembly}}.bam", **config)
    output:
        bedgraphish=expand("{result_dir}/genrich/{{sample}}-{{assembly}}.bdgish", **config),
        log=expand("{result_dir}/genrich/{{sample}}-{{assembly}}.log", **config)
    conda:
        "../envs/call_peak.yaml"
    params:
        config['peak_caller'].get('genrich', " ")
    threads: 15
    shell:
        "input=$(echo {input} | tr ' ' ','); "
        "Genrich -X -t $input -f {output.log} -k {output.bedgraphish} {params}"


rule call_peak_genrich:
    # TODO log
    input:
        log=expand("{result_dir}/genrich/{{sample}}-{{assembly}}.log", **config)
    output:
        narrowpeak= expand("{result_dir}/genrich/{{sample}}-{{assembly}}_peaks.narrowPeak", **config)
    conda:
        "../envs/call_peak.yaml"
    params:
        config['peak_caller']['genrich']
    threads: 1
    shell:
        "Genrich -P -f {input.log} -o {output.narrowpeak} {params}"


config['macs2_types'] = ['control_lambda.bdg', 'summits.bed', 'peaks.narrowPeak',
                         'peaks.xls', 'treat_pileup.bdg']
rule call_peak_macs2:
    #
    # Calculates genome size based on unique kmers of average length
    #
    input:
        bam=   expand("{result_dir}/dedup/{{sample}}-{{assembly}}.bam", **config),
        #fastqc=expand("{result_dir}/trimmed/{{sample}}-{{condition}}-{{project}}-{{assembly}}_trimmed_fastqc.zip", **config)
    output:
        expand("{result_dir}/macs2/{{sample}}-{{assembly}}_{macs2_types}", **config)
    params:
        name="{sample}-{assembly}",
        genome=f"{config['genome_dir']}/{{assembly}}/{{assembly}}.fa",
        macs_params=config['peak_caller']['macs2']
    log:
        "logs/call_peak_macs2/{sample}-{assembly}.log"
    conda:
        "../envs/call_peak_macs2.yaml"
    shell:
# TODO: calculate genome size
#         kmer_size=$(unzip -p {{input.fastqc}} {{params.name}}_trimmed_fastqc/fastqc_data.txt  | grep -P -o '(?<=Sequence length\\t).*' | grep -P -o '\d+$');
#         GENSIZE=$(unique-kmers.py {{params.genome}} -k $kmer_size --quiet 2>&1 | grep -P -o '(?<=\.fa: ).*');
        f"""
        GENSIZE=123456789
        macs2 callpeak -t {{input.bam}} --outdir {config['result_dir']}/macs2/ -n {{params.name}} \
        {{params.macs_params}} -g $GENSIZE > {{log}} 2>&1
        """
