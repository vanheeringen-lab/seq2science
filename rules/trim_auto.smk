# ruleorder in favour of PE over SE, otherwise can't resolve properly
ruleorder: trim_galore_PE > trim_galore_SE

rule trim_galore_SE:
    input:
        download_fastq(expand("{result_dir}/{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config))
    output:
        se=expand("{result_dir}/{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config),
        qc=expand("{result_dir}/{trimmed_dir}/{{sample}}.{fqsuffix}.gz_trimming_report.txt", **config)
    conda:
        "../envs/trim_auto.yaml"
    threads: 6
    log:
        expand("{log_dir}/trim_galore/{{sample}}.log", **config)  
    params:
        config=config['trim_galore'],
        fqsuffix=config['fqsuffix']
    shell:
        """
        cpulimit --include-children -l {threads}00 --\
        trim_galore -j {threads} {params.config} -o $(dirname {output.se}) {input} > {log} 2>&1

        # now rename to proper output
        rename 's/.fq/.{params.fqsuffix}/' "$(dirname {output.se})/{wildcards.sample}_trimmed.fq.gz"
        """


rule trim_galore_PE:
    input:
        r1=download_fastq(expand("{result_dir}/{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config)),
        r2=download_fastq(expand("{result_dir}/{fastq_dir}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config))
    output:
        r1=expand("{result_dir}/{trimmed_dir}/{{sample}}_{fqext1}_trimmed.{fqsuffix}.gz", **config),
        r2=expand("{result_dir}/{trimmed_dir}/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config),
        qc=expand("{result_dir}/{trimmed_dir}/{{sample}}_{fqext}.{fqsuffix}.gz_trimming_report.txt", **config)
    conda:
        "../envs/trim_auto.yaml"
    threads: 6
    log:
        expand("{log_dir}/trim_galore/{{sample}}.log", **config)
    params:
        config=config['trim_galore'],
        fqsuffix=config['fqsuffix']
    shell:
        """
        cpulimit --include-children -l {threads}00 --\
        trim_galore --paired -j {threads} {params.config} -o $(dirname {output.r1}) {input.r1} {input.r2} > {log} 2>&1

        # now rename to proper output
        find "$(dirname {output.r1})/" -name "{wildcards.sample}_*val_*.fq.gz" | rename 's/.fq/.{params.fqsuffix}/'
        find "$(dirname {output.r1})/" -name "{wildcards.sample}_*val_*.fastq.gz" | rename 's/_val_\d/_trimmed/'
        """
