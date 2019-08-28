ruleorder: trim_galore_PE > trim_galore_SE

rule trim_galore_SE:
    """
    Automated adapter detection, adapter trimming, and quality trimming through trim galore (single-end).
    """
    input:
        expand("{result_dir}/{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config)
    output:
        se=expand("{result_dir}/{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config),
        qc=expand("{result_dir}/{qc_dir}/trimming/{{sample}}.{fqsuffix}.gz_trimming_report.txt", **config)
    conda:
        "../envs/trimgalore.yaml"
    threads: 6
    log:
        expand("{log_dir}/trim_galore_SE/{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/trim_galore_SE/{{sample}}.benchmark.txt", **config)[0]
    params:
        config=config['trim_galore'],
        fqsuffix=config['fqsuffix']
    wildcard_constraints:
        sample=f".*(?<!{config['fqext1']}|{config['fqext2']})+"
    shell:
        """
        cpulimit --include-children -l {threads}00 --\
        trim_galore -j {threads} {params.config} -o $(dirname {output.se}) {input} > {log} 2>&1

        # now rename to proper output
        if [ "{params.fqsuffix}" != "fq" ]; then
          mv "$(dirname {output.se})/{wildcards.sample}_trimmed.fq.gz" {output.se}
        fi 
        
        # move the trimming report to the desired directory
        report=$(dirname {output.se})/{wildcards.sample}.{params.fqsuffix}.gz_trimming_report.txt
        if [[ -f $report ]]; then
            mv $report {output.qc}
        fi
        """


rule trim_galore_PE:
    """
    Automated adapter detection, adapter trimming, and quality trimming through trim galore (paired-end).
    """
    input:
        r1=expand("{result_dir}/{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config),
        r2=expand("{result_dir}/{fastq_dir}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config)
    output:
        r1=expand("{result_dir}/{trimmed_dir}/{{sample}}_{fqext1}_trimmed.{fqsuffix}.gz", **config),
        r2=expand("{result_dir}/{trimmed_dir}/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config),
        qc=expand("{result_dir}/{qc_dir}/trimming/{{sample}}_{fqext}.{fqsuffix}.gz_trimming_report.txt", **config)
    conda:
        "../envs/trimgalore.yaml"
    threads: 6
    log:
        expand("{log_dir}/trim_galore_PE/{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/trim_galore_PE/{{sample}}.benchmark.txt", **config)[0]
    params:
        config=config['trim_galore'],
        fqsuffix=config['fqsuffix'] # ,
        # fqext=config['fqext']
    shell:
        """
        cpulimit --include-children -l {threads}00 --\
        trim_galore --paired -j {threads} {params.config} -o $(dirname {output.r1}) {input.r1} {input.r2} > {log} 2>&1

        # now rename to proper output
        for f in $(find "$(dirname {output.r1})/" -name "{wildcards.sample}_*val_*.fq.gz"); do
            mv "$f" "$(echo "$f" | sed s/.fq/.{params.fqsuffix}/)"; done
        for f in $(find "$(dirname {output.r1})/" -maxdepth 1 -name "{wildcards.sample}_*val_*.{params.fqsuffix}.gz"); do
            mv "$f" "$(echo "$f" | sed s/_val_./_trimmed/)"; done

        # move the trimming reports
        for f in $(find "$(dirname {output.r1})/" -name "{wildcards.sample}_*.{params.fqsuffix}.gz_trimming_report.txt"); do
            mv "$f" "$(dirname {output.qc[0]})/$(basename $f)"; done
        """

        # touch $(dirname {output.r1})/{wildcards.sample}_pass_1.fastq.gz_trimming_report.txt
        # touch $(dirname {output.r1})/{wildcards.sample}_pass_2.fastq.gz_trimming_report.txt
        # for f in $(find "$(dirname {output.r1})/" -name "{wildcards.sample}_*.{params.fqsuffix}.gz_trimming_report.txt"); do
        #     mv "$f" "$(dirname {output.qc[0]})/$(basename $f)"; done
        #     echo $f
        #     echo $(dirname {output.qc[0]})/$(basename $f)
        #     echo ''
        # done
        # exit

        # report1=$(dirname {output.r1})/{wildcards.sample}_{params.fqext[0]}.{params.fqsuffix}.gz_trimming_report.txt
        # report2=$(dirname {output.r1})/{wildcards.sample}_{params.fqext[1]}.{params.fqsuffix}.gz_trimming_report.txt
        # if [[ -f $report1 ]]; then
        #     mv $report1 {output.qc[0]}
        #     mv $report2 {output.qc[1]}
        # fi

if 'condition' in samples and config.get('combine_replicates', '') == 'merge':
    def get_merge_replicates(wildcards):
        return expand([f"{{result_dir}}/{{trimmed_dir}}/{replicate}{wildcards.fqext}_trimmed.{{fqsuffix}}.gz"
               for replicate in samples[samples['condition'] == wildcards.condition].index], **config)

    rule merge_replicates:
        """
        Merge replicates (fastqs) simply by concatenating the files.
        """
        input:
            get_merge_replicates
        output:
            sorted(expand("{result_dir}/{trimmed_dir}/merged/{{condition}}{{fqext}}_trimmed.{fqsuffix}.gz", **config))
        wildcard_constraints:
            fqext=".*",
            condition="[^/_]*"
        log:
            expand("{log_dir}/merge_replicates/{{condition}}{{fqext}}.log", **config)
        benchmark:
            expand("{benchmark_dir}/merge_replicates/{{condition}}{{fqext}}.benchmark.txt", **config)[0]
        shell:
            "cat {input} > {output} 2> {log}"
