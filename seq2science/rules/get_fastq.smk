import glob
import os
import re


rule run2sra:
    """
    Download the SRA of a sample by its unique identifier.

    Tries first downloading with the faster ascp protocol, if that fails it
    falls back on the slower http protocol.
    """
    output:
        temp(directory(expand("{sra_dir}/{{run}}", **config))),
    log:
        expand("{log_dir}/run2sra/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/run2sra/{{run}}.benchmark.txt", **config)[0]
    message: explain_rule("run2sra")
    resources:
        parallel_downloads=1,
    conda:
        "../envs/get_fastq.yaml"
    wildcard_constraints:
        run="SRR\d+",
    shell:
        """
        # three attempts
        for i in {{1..3}}
        do
            # acquire a lock
            (
                flock --timeout 30 200 || continue
                sleep 2
            ) 200>{eutils_cache_lock}

            # dump
            prefetch --max-size 999999999999 --output-directory {output} --log-level debug --progress {wildcards.run} >> {log} 2>&1 && break
            sleep 10
        done
        """


rule sra2fastq_SE:
    """
    Downloaded (raw) SRAs are converted to single-end fastq files.
    """
    input:
        rules.run2sra.output,
    output:
        fastq=temp(expand("{fastq_dir}/runs/{{run}}.{fqsuffix}.gz", **config)),
        tmpdir=temp(directory(expand("{fastq_dir}/runs/tmp/{{run}}", **config))),
    log:
        expand("{log_dir}/sra2fastq_SE/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sra2fastq_SE/{{run}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        run="SRR\d+",
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        """
        # acquire the lock
        (
            flock --timeout 30 200 || exit 1
            sleep 3
        ) 200>{eutils_cache_lock}

        # setup tmp dir
        mkdir -p {output.tmpdir}

        # dump to tmp dir
        parallel-fastq-dump -s {input}/* -O {output.tmpdir} \
        --threads {threads} --split-spot --skip-technical --dumpbase --readids \
        --clip --read-filter pass --defline-seq '@$ac.$si.$sg/$ri' \
        --defline-qual '+' --gzip >> {log} 2>&1

        # rename file and move to output dir
        mv {output.tmpdir}/* {output.fastq}
        """


rule sra2fastq_PE:
    """
    Downloaded (raw) SRAs are converted to paired-end fastq files.
    Forward and reverse samples will be switched if forward/reverse names are not lexicographically ordered.
    """
    input:
        rules.run2sra.output,
    output:
        fastq=temp(expand("{fastq_dir}/runs/{{run}}_{fqext}.{fqsuffix}.gz", **config)),
        tmpdir=temp(directory(expand("{fastq_dir}/runs/tmp/{{run}}", **config))),
    log:
        expand("{log_dir}/sra2fastq_PE/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sra2fastq_PE/{{run}}.benchmark.txt", **config)[0]
    threads: 8
    wildcard_constraints:
        run="SRR\d+",
    conda:
        "../envs/get_fastq.yaml"
    shell:
        """
        # acquire the lock
        (
            flock --timeout 30 200 || exit 1
            sleep 3
        ) 200>{eutils_cache_lock}

        # setup tmp dir
        mkdir -p {output.tmpdir}

        # dump to tmp dir
        parallel-fastq-dump -s {input}/* -O {output.tmpdir} \
        --threads {threads} --split-e --skip-technical --dumpbase \
        --readids --clip --read-filter pass --defline-seq '@$ac.$si.$sg/$ri' \
        --defline-qual '+' --gzip >> {log} 2>&1

        # rename file and move to output dir
        mv {output.tmpdir}/*_1* {output.fastq[0]}
        mv {output.tmpdir}/*_2* {output.fastq[1]}
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
        temp(expand("{fastq_dir}/runs/{{run}}.{fqsuffix}.gz", **config)),
    resources:
        parallel_downloads=1,
    log:
        expand("{log_dir}/ena2fastq_SE/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/ena2fastq_SE/{{run}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        run="|".join(ena_single_end) if len(ena_single_end) else "$a"
    run:
        shell("mkdir -p {config[fastq_dir]}/tmp/ >> {log} 2>&1")
        url = run2download[wildcards.run]
        if config.get('ascp_path') and config.get('ascp_key'):
            shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {url} {output} >> {log} 2>&1")
        else:
            shell("wget {url} -O {output} >> {log} 2>&1")


rule ena2fastq_PE:
    """
    Download paired-end fastq files directly from the ENA.
    """
    output:
        temp(expand("{fastq_dir}/runs/{{run}}_{fqext}.{fqsuffix}.gz", **config)),
    resources:
        parallel_downloads=1,
    log:
        expand("{log_dir}/ena2fastq_PE/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/ena2fastq_PE/{{run}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        run="|".join(ena_paired_end) if len(ena_paired_end) else "$a"
    run:
        shell("mkdir -p {config[fastq_dir]}/tmp >> {log} 2>&1")
        urls = run2download[wildcards.run]
        if config.get('ascp_path') and config.get('ascp_key'):
            shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {urls[0]} {output[0]} >> {log} 2>&1")
            shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {urls[1]} {output[1]} >> {log} 2>&1")
        else:
            shell("wget {urls[0]} -O {output[0]} >> {log} 2>&1")
            shell("wget {urls[1]} -O {output[1]} >> {log} 2>&1")


def get_runs_from_sample(wildcards):
    run_fastqs = []
    for run in sampledict[wildcards.sample].get("runs", []):
        run_fastqs.append(f"{config['fastq_dir']}/runs/{run}{wildcards.suffix}")

    return run_fastqs


rule runs2sample:
    """
    Concatenate a single run or multiple runs together into a fastq
    """
    input:
        get_runs_from_sample
    output:
        expand("{fastq_dir}/{{sample}}{{suffix}}", **config),
    log:
        expand("{log_dir}/run2sample/{{sample}}{{suffix}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/run2sample/{{sample}}{{suffix}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        run="(GSM|SRR|ERR|DRR)\d+",
    run:
        shell("cp {input[0]} {output}")

        # now append all the later ones
        for i in range(1, len(input)):
            inputfile = input[i]
            shell("cat {inputfile} >> {output}")
