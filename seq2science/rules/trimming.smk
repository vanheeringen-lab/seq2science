if config["trimmer"] == "trimgalore":
    ruleorder: trimgalore_PE > trimgalore_SE


    rule trimgalore_SE:
        """
        Automated adapter detection, adapter trimming, and quality trimming through trim galore (single-end).
        """
        input:
            expand("{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config),
        output:
            se=temp(expand("{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config)),
            qc=expand("{qc_dir}/trimming/{{sample}}.{fqsuffix}.gz_trimming_report.txt", **config),
        conda:
            "../envs/trimgalore.yaml"
        threads: 6
        message: explain_rule("trim_galore_SE")
        log:
            expand("{log_dir}/trim_galore_SE/{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/trim_galore_SE/{{sample}}.benchmark.txt", **config)[0]
        params:
            config=config["trimoptions"],
            fqsuffix=config["fqsuffix"],
        shell:
            ("cpulimit --include-children -l {threads}00 --" if config. get("cpulimit", True) else"")+
            """\
            trim_galore -j {threads} {params.config} -o $(dirname {output.se}) {input} > {log} 2>&1
    
            # now rename to proper output
            if [ "{params.fqsuffix}" != "fq" ]; then
              mv "$(dirname {output.se})/{wildcards.sample}_trimmed.fq.gz" {output.se}
            fi 
    
            # move the trimming report to qc directory
            report=$(dirname {output.se})/{wildcards.sample}.{params.fqsuffix}.gz_trimming_report.txt
            mv $report {output.qc}
            """


    rule trimgalore_PE:
        """
        Automated adapter detection, adapter trimming, and quality trimming through trim galore (paired-end).
        """
        input:
            r1=expand("{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config),
            r2=expand("{fastq_dir}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config),
        output:
            r1=temp(expand("{trimmed_dir}/{{sample}}_{fqext1}_trimmed.{fqsuffix}.gz", **config)),
            r2=temp(expand("{trimmed_dir}/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config)),
            qc=expand("{qc_dir}/trimming/{{sample}}_{fqext}.{fqsuffix}.gz_trimming_report.txt", **config),
        conda:
            "../envs/trimgalore.yaml"
        threads: 6
        message: explain_rule("trimgalore_PE")
        log:
            expand("{log_dir}/trimgalore_PE/{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/trimgalore_PE/{{sample}}.benchmark.txt", **config)[0]
        params:
            config=config["trimoptions"],
            fqsuffix=config["fqsuffix"],
        shell:
            ("cpulimit --include-children -l {threads}00 --" if config. get("cpulimit", True) else"")+
            """\
            trim_galore --paired -j {threads} {params.config} -o $(dirname {output.r1}) {input.r1} {input.r2} > {log} 2>&1
    
            # now rename to proper output
            for f in $(find "$(dirname {output.r1})/" -name "{wildcards.sample}_*val_*.fq.gz"); do
                mv "$f" "$(echo "$f" | sed s/.fq/.{params.fqsuffix}/)"; done
            for f in $(find "$(dirname {output.r1})/" -maxdepth 1 -name "{wildcards.sample}_*val_*.{params.fqsuffix}.gz"); do
                mv "$f" "$(echo "$f" | sed s/_val_./_trimmed/)"; done
    
            # move the trimming reports to qc directory
            for f in $(find "$(dirname {output.r1})/" -name "{wildcards.sample}_*.{params.fqsuffix}.gz_trimming_report.txt"); do
                mv "$f" "$(dirname {output.qc[0]})/$(basename $f)"; done
            """


elif config["trimmer"] == "fastp":
    ruleorder: fastp_PE > fastp_SE


    rule fastp_SE:
        """
        Automated adapter detection, adapter trimming, and quality trimming through fastp (single-end).
        """
        input:
            expand("{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config),
        output:
            se=temp(expand("{trimmed_dir}/{{sample}}_trimmed.{fqsuffix}.gz", **config)),
            qc_json=expand("{qc_dir}/trimming/{{sample}}.fastp.json", **config),
            qc_html=expand("{qc_dir}/trimming/{{sample}}.fastp.html", **config),
        conda:
            "../envs/fastp.yaml"
        threads: 3
        message: explain_rule("fastp_SE")
        wildcard_constraints:
            sample="|".join([sample for sample in samples.index if config["layout"][sample] == "SINGLE"])
        log:
            expand("{log_dir}/fastp_SE/{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/fastp_SE/{{sample}}.benchmark.txt", **config)[0]
        params:
            fqsuffix=config["fqsuffix"],
        shell:
            """\
            fastp -w {threads} --in1 {input} --out1 {output.se} -h {output.qc_html} -j {output.qc_json} \
            > {log} 2>&1
            """


    rule fastp_PE:
        """
        Automated adapter detection, adapter trimming, and quality trimming through fastp (paired-end).
        """
        input:
            r1=expand("{fastq_dir}/{{sample}}_{fqext1}.{fqsuffix}.gz", **config),
            r2=expand("{fastq_dir}/{{sample}}_{fqext2}.{fqsuffix}.gz", **config),
        output:
            r1=temp(expand("{trimmed_dir}/{{sample}}_{fqext1}_trimmed.{fqsuffix}.gz", **config)),
            r2=temp(expand("{trimmed_dir}/{{sample}}_{fqext2}_trimmed.{fqsuffix}.gz", **config)),
            qc_json=expand("{qc_dir}/trimming/{{sample}}.fastp.json", **config),
            qc_html=expand("{qc_dir}/trimming/{{sample}}.fastp.html", **config),
        conda:
            "../envs/fastp.yaml"
        threads: 3
        wildcard_constraints:
            sample="|".join([sample for sample in samples.index if config["layout"][sample] == "PAIRED"])
        message: explain_rule("fastp_PE")
        log:
            expand("{log_dir}/fastp_PE/{{sample}}.log", **config),
        benchmark:
            expand("{benchmark_dir}/fastp_PE/{{sample}}.benchmark.txt", **config)[0]
        params:
            config=config["trimoptions"],
        shell:
            """\
            fastp -w {threads} --in1 {input[0]} --in2 {input[1]} \
            --out1 {output.r1} --out2 {output.r2} -h {output.qc_html} -j {output.qc_json} \
            {params.config} > {log} 2>&1
            """
