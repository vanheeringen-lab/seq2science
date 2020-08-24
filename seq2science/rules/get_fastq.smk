import glob
import os
import re


def gsm2srx(gsm):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"

    conn = urllib.request.urlopen(url)
    html = conn.read()

    soup = BeautifulSoup(html, features="html.parser")
    links = soup.find_all('a')

    for tag in links:
        link = tag.get('href', None)
        if link is not None and 'SRX' in link:
            SRX = link[link.find("SRX"):]
            return SRX
    raise ValueError(f"Sample {gsm} has been put in wrongly in the SRA and "
                     f"seq2science is not capable of downloading it...")


rule id2sra:
    """
    Download the SRA of a sample by its unique identifier.

    Tries first downloading with the faster ascp protocol, if that fails it 
    falls back on the slower http protocol.
    """
    output:
        temp(directory(expand("{sra_dir}/{{sample}}", **config))),
    log:
        expand("{log_dir}/id2sra/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/id2sra/{{sample}}.benchmark.txt", **config)[0]
    message: explain_rule("id2sra")
    resources:
        parallel_downloads=1,
    conda:
        "../envs/get_fastq.yaml"
    wildcard_constraints:
        sample="(GSM|SRR|ERR|DRR)\d+",
    shell:
        """
        # three attempts
        for i in {{1..3}}
        do
            # acquire a lock
            (
                flock --timeout 30 200 || continue
                sleep 2
            ) 200>{layout_cachefile_lock}
    
            # dump
            prefetch --max-size 999999999999 --output-directory {output} --log-level debug --progress {wildcards.sample} >> {log} 2>&1 && break
            sleep 10
        done
        """



rule sra2fastq_SE:
    """
    Downloaded (raw) SRAs are converted to single-end fastq files.
    """
    input:
        rules.id2sra.output,
    output:
        fastq=expand("{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config),
        tmp_dump=temp(directory(expand("{sra_dir}/tmp/{{sample}}", **config))),
        tmp_fastq=temp(directory(expand("{sra_dir}/fastq/{{sample}}", **config))),
    log:
        expand("{log_dir}/sra2fastq_SE/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sra2fastq_SE/{{sample}}.benchmark.txt", **config)[0]
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        """
        # acquire a lock
        (
            flock --timeout 30 200 || exit 1
            sleep 2
        ) 200>{layout_cachefile_lock}

        # dump
        fasterq-dump -s {input}/* -O {output.tmp_fastq} -t {output.tmp_dump} --threads {threads} --split-spot >> {log} 2>&1

        # rename file and move to output dir
        for f in $(ls -1q {output.tmp_fastq} | grep -oP "^[^_]+" | uniq); do
            dst={config[fastq_dir]}/{wildcards.sample}.{config[fqsuffix]}
            cat "{output.tmp_fastq}/${{f}}" >> $dst
        done
        pigz -p {threads} {config[fastq_dir]}/{wildcards.sample}.{config[fqsuffix]}
        """


rule sra2fastq_PE:
    """
    Downloaded (raw) SRAs are converted to paired-end fastq files.
    Forward and reverse samples will be switched if forward/reverse names are not lexicographically ordered.
    """
    input:
        rules.id2sra.output,
    output:
        fastq=expand("{fastq_dir}/{{sample}}_{fqext}.{fqsuffix}.gz", **config),
        tmp_dump=temp(directory(expand("{sra_dir}/tmp/{{sample}}", **config))),
        tmp_fastq=temp(directory(expand("{sra_dir}/fastq/{{sample}}", **config))),
    log:
        expand("{log_dir}/sra2fastq_PE/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sra2fastq_PE/{{sample}}.benchmark.txt", **config)[0]
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        """
        # acquire the lock
        (
            flock --timeout 30 200 || exit 1
            sleep 2
        ) 200>{layout_cachefile_lock}

        # dump
        fasterq-dump -s {input}/* -O {output.tmp_fastq} -t {output.tmp_dump} --threads {threads} --split-3 >> {log} 2>&1

        # rename files and move to output dir
        for f in $(ls -1q {output.tmp_fastq} | grep -oP "^[^_]+" | uniq); do
            dst_1={config[fastq_dir]}/{wildcards.sample}_{config[fqext1]}.{config[fqsuffix]}
            dst_2={config[fastq_dir]}/{wildcards.sample}_{config[fqext2]}.{config[fqsuffix]}
            cat "{output.tmp_fastq}/${{f}}_1.fastq" >> $dst_1
            cat "{output.tmp_fastq}/${{f}}_2.fastq" >> $dst_2
        done
        pigz -p {threads} {config[fastq_dir]}/{wildcards.sample}_{config[fqext1]}.{config[fqsuffix]}
        pigz -p {threads} {config[fastq_dir]}/{wildcards.sample}_{config[fqext2]}.{config[fqsuffix]}
        """


# NOTE: if the workflow fails it tends to blame this rule.
# Set "debug: True" in the config to see the root cause.
if not config.get("debug"):
    ruleorder: renamefastq_PE > sra2fastq_PE

    def get_wrong_fqext(wildcards):
        """get all local samples with fqexts that do not match the config"""
        fastqs = glob.glob(os.path.join(config["fastq_dir"], wildcards.sample + "*" + config["fqsuffix"] + ".gz"))
        # exclude samples with the correct fqext
        r1 = wildcards.sample + "_" + config["fqext1"] + "." + config["fqsuffix"] + ".gz"
        r2 = wildcards.sample + "_" + config["fqext2"] + "." + config["fqsuffix"] + ".gz"
        fastqs = [f for f in fastqs if not re.match(r1 + "|" + r2, f)]
        if len(fastqs) == 0:
            fastqs.append(
                f"If you can read this, seq2science is looking for non-existing files (for sample {wildcards}) or in the wrong location. "
                "Tip: Try removing the seq2science cache with 'seq2science clean'"
            )
        return sorted(fastqs)


    rule renamefastq_PE:
        """
        Create symlinks to fastqs with incorrect fqexts (default R1/R2).
        Forward and reverse samples will be switched if forward/reverse names are not 
        lexicographically ordered.
        """
        input:
            get_wrong_fqext,
        output:
            temp(expand("{fastq_dir}/{{sample}}_{fqext}.{fqsuffix}.gz", **config)),
        shell:
            """
            ln {input[0]} {output[0]}
            ln {input[1]} {output[1]}
            """


ruleorder: ena2fastq_SE > sra2fastq_SE
ruleorder: ena2fastq_PE > sra2fastq_PE
ruleorder: ena2fastq_PE > ena2fastq_SE


rule ena2fastq_SE:
    """
    Download single-end fastq files directly from the ENA.
    """
    output:
        expand("{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config),
    resources:
        parallel_downloads=1,
    log:
        expand("{log_dir}/ena2fastq_SE/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/ena2fastq_SE/{{sample}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        sample="|".join(ena_single_end_urls.keys()) if len(ena_single_end_urls) else "$a"
    run:
        try:
            shell("mkdir -p {config[fastq_dir]}/{wildcards.sample} >> {log} 2>&1")
            for srr, url in ena_single_end_urls[wildcards.sample]:
                if config.get('ascp_path') and config.get('ascp_key'):
                    shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {url} {config[fastq_dir]}/{wildcards.sample}/{srr}.{config[fqsuffix]}.gz >> {log} 2>&1")
                else:
                    shell("wget {url} -O {config[fastq_dir]}/{wildcards.sample}/{srr}.{config[fqsuffix]}.gz >> {log} 2>&1")
                shell("cat {config[fastq_dir]}/{wildcards.sample}/{srr}.{config[fqsuffix]}.gz >> {output} 2> {log}")
        except:
            pass
        finally:
            shell("rm -r {config[fastq_dir]}/{wildcards.sample}/ >> {log} 2>&1")


rule ena2fastq_PE:
    """
    Download paired-end fastq files directly from the ENA.
    """
    output:
        expand("{fastq_dir}/{{sample}}_{fqext}.{fqsuffix}.gz", **config),
    resources:
        parallel_downloads=1,
    log:
        expand("{log_dir}/ena2fastq_PE/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/ena2fastq_PE/{{sample}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        sample="|".join(ena_paired_end_urls.keys()) if len(ena_paired_end_urls) else "$a"
    run:
        try:
            shell("mkdir -p {config[fastq_dir]}/{wildcards.sample}/ >> {log} 2>&1")
            for srr, urls in ena_paired_end_urls[wildcards.sample]:
                if config.get('ascp_path') and config.get('ascp_key'):
                    shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {urls[0]} {config[fastq_dir]}/{wildcards.sample}/{srr}_{config[fqext1]}.{config[fqsuffix]}.gz >> {log} 2>&1")
                    shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {urls[1]} {config[fastq_dir]}/{wildcards.sample}/{srr}_{config[fqext2]}.{config[fqsuffix]}.gz >> {log} 2>&1")
                else:
                    shell("wget {urls[0]} -O {config[fastq_dir]}/{wildcards.sample}/{srr}_{config[fqext1]}.{config[fqsuffix]}.gz >> {log} 2>&1")
                    shell("wget {urls[1]} -O {config[fastq_dir]}/{wildcards.sample}/{srr}_{config[fqext2]}.{config[fqsuffix]}.gz >> {log} 2>&1")
                shell("cat {config[fastq_dir]}/{wildcards.sample}/{srr}_{config[fqext1]}.{config[fqsuffix]}.gz >> {output[0]} 2> {log}")
                shell("cat {config[fastq_dir]}/{wildcards.sample}/{srr}_{config[fqext2]}.{config[fqsuffix]}.gz >> {output[1]} 2> {log}")
        except:
            pass
        finally:
            shell("rm -r {config[fastq_dir]}/{wildcards.sample}/ >> {log} 2>&1")
