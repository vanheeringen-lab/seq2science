"""
all rules/logic related downloading public fastqs should be here.
"""
localrules: ena2fastq_SE, ena2fastq_PE, gsa2fastq_SE, gsa2fastq_PE, runs2sample


import os


# ENA > SRA
ruleorder: ena2fastq_SE > gsa2fastq_SE > sra2fastq_SE
ruleorder: ena2fastq_PE > gsa2fastq_PE > sra2fastq_PE
# PE > SE
ruleorder: ena2fastq_PE > ena2fastq_SE
ruleorder: sra2fastq_PE > sra2fastq_SE
ruleorder: gsa2fastq_PE > gsa2fastq_SE


rule sra2fastq_SE:
    """
    Downloaded (raw) SRAs are converted to single-end fastq files.
    """
    output:
        temp(expand("{fastq_dir}/runs/{{run}}.{fqsuffix}.gz", **config))
    log:
        expand("{log_dir}/sra2fastq_SE/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sra2fastq_SE/{{run}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        run="|".join(SRA_SINGLE_END) if len(SRA_SINGLE_END) else "$a",  # only try to dump (single-end) SRA samples
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    retries: 2
    shell:
        """
        tmpdir=`mktemp -d`
        echo "Using tmpdir $tmpdir" >> {log}

        # move to output dir since somehow prefetch and parallel-fastq-dump sometimes put files in the cwd...
        cd $tmpdir

        # acquire the lock
        (
            flock -w 30 200 || exit 1
            sleep 3
        ) 200>{PYSRADB_CACHE_LOCK}

        # three attempts
        for i in {{1..3}}
        do
            # dump
            source {workflow.basedir}/../../scripts/download_sra.sh {wildcards.run} {log} && break
            sleep 10
        done

        # dump to tmp dir
        parallel-fastq-dump -s {wildcards.run}/{wildcards.run}.sra -O $tmpdir \
        --threads {threads} --split-spot --skip-technical --dumpbase --readids \
        --clip --read-filter pass --defline-seq '@$ac.$si.$sg/$ri' \
        --defline-qual '+' --gzip >> {log} 2>&1

        # rename file and move to output dir
        mv $tmpdir/*fastq.gz {output} >> {log} 2>&1

        # Remove temporary directory
        rm -rf $tmpdir >> {log} 2>&1
        """


rule sra2fastq_PE:
    """
    Downloaded (raw) SRAs are converted to paired-end fastq files.
    Forward and reverse samples will be switched if forward/reverse names are not lexicographically ordered.
    """
    output:
        temp(expand("{fastq_dir}/runs/{{run}}_{fqext}.{fqsuffix}.gz", **config))
    log:
        expand("{log_dir}/sra2fastq_PE/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sra2fastq_PE/{{run}}.benchmark.txt", **config)[0]
    threads: 8
    wildcard_constraints:
        run="|".join(SRA_PAIRED_END) if len(SRA_PAIRED_END) else "$a",  # only try to dump (paired-end) SRA samples
    conda:
        "../envs/get_fastq.yaml"
    retries: 2
    shell:
        """
        tmpdir=`mktemp -d`
        echo "Using tmpdir $tmpdir" >> {log}

        # move to output dir since somehow prefetch and parallel-fastq-dump sometimes put files in the cwd...
        cd $tmpdir

        # acquire the lock
        (
            flock -w 30 200 || exit 1
            sleep 3
        ) 200>{PYSRADB_CACHE_LOCK}

        # three attempts
        for i in {{1..3}}
        do
            # dump
            source {workflow.basedir}/../../scripts/download_sra.sh {wildcards.run} {log} && break
            sleep 10
        done

        # rename file and move to output dir
        mv $tmpdir/*_1* {output.fastq[0]} >> {log} 2>&1
        mv $tmpdir/*_2* {output.fastq[1]} >> {log} 2>&1
        
        # Remove temporary directory
        rm -rf $tmpdir >> {log} 2>&1
        """

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
        run="|".join(ENA_SINGLE_END) if len(ENA_SINGLE_END) else "$a",
    retries: 2
    run:
        shell("mkdir -p {config[fastq_dir]}/tmp/ >> {log} 2>&1")
        url = RUN2DOWNLOAD[wildcards.run]
        if config.get("ascp_path") and config.get("ascp_key"):
            shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {url} {output} >> {log} 2>&1")
        else:
            shell("wget {url} -O {output} --waitretry 20 >> {log} 2>&1")

        if not (os.path.exists(output[0]) and os.path.getsize(output[0]) > 0):
            shell('echo "Something went wrong.. The downloaded file was empty!" >> {log} 2>&1')
            shell("exit 1 >> {log} 2>&1")


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
        run="|".join(ENA_PAIRED_END) if len(ENA_PAIRED_END) else "$a",
    retries: 2
    run:
        shell("mkdir -p {config[fastq_dir]}/tmp >> {log} 2>&1")
        urls = RUN2DOWNLOAD[wildcards.run]
        if config.get("ascp_path") and config.get("ascp_key"):
            shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {urls[0]} {output[0]} >> {log} 2>&1")
            shell("{config[ascp_path]} -QT -l 1G -P33001 -i {config[ascp_key]} {urls[1]} {output[1]} >> {log} 2>&1")
        else:
            shell("wget {urls[0]} -O {output[0]} --waitretry 20 >> {log} 2>&1")
            shell("wget {urls[1]} -O {output[1]} --waitretry 20 >> {log} 2>&1")

        if not (
            os.path.exists(output[0])
            and (os.path.getsize(output[0]) > 0)
            and os.path.exists(output[1])
            and (os.path.getsize(output[1]) > 0)
        ):
            shell('echo "Something went wrong.. The downloaded file(s) were empty!" >> {log} 2>&1')
            shell("exit 1 >> {log} 2>&1")


rule gsa2fastq_SE:
    """
    Download single-end fastq files from the GSA.
    """
    output:
        temp(expand("{fastq_dir}/runs/{{run}}.{fqsuffix}.gz", **config)),
    resources:
        parallel_downloads=1,
    log:
        expand("{log_dir}/gsa2fastq_SE/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/gsa2fastq_SE/{{run}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        run="|".join(GSA_SINGLE_END) if len(GSA_SINGLE_END) else "$a",
    retries: 2
    run:
        shell("mkdir -p {config[fastq_dir]}/tmp/ >> {log} 2>&1")
        url = RUN2DOWNLOAD[wildcards.run]

        shell("wget {url} -O {output} --waitretry 20 >> {log} 2>&1")

        if not (os.path.exists(output[0]) and os.path.getsize(output[0]) > 0):
            shell('echo "Something went wrong.. The downloaded file was empty!" >> {log} 2>&1')
            shell("exit 1 >> {log} 2>&1")


rule gsa2fastq_PE:
    """
    Download paired-end fastq files from the GSA.
    """
    output:
        temp(expand("{fastq_dir}/runs/{{run}}_{fqext}.{fqsuffix}.gz", **config)),
    resources:
        parallel_downloads=1,
    log:
        expand("{log_dir}/gsa2fastq_PE/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/gsa2fastq_PE/{{run}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        run="|".join(GSA_PAIRED_END) if len(GSA_PAIRED_END) else "$a",
    retries: 2
    run:
        shell("mkdir -p {config[fastq_dir]}/tmp >> {log} 2>&1")
        urls = RUN2DOWNLOAD[wildcards.run]

        shell("wget {urls[0]} -O {output[0]} --waitretry 20 >> {log} 2>&1")
        shell("wget {urls[1]} -O {output[1]} --waitretry 20 >> {log} 2>&1")

        if not (
            os.path.exists(output[0])
            and (os.path.getsize(output[0]) > 0)
            and os.path.exists(output[1])
            and (os.path.getsize(output[1]) > 0)
        ):
            shell('echo "Something went wrong.. The downloaded file(s) were empty!" >> {log} 2>&1')
            shell("exit 1 >> {log} 2>&1")


def get_runs_from_sample(wildcards):
    run_fastqs = []
    for run in SAMPLEDICT[wildcards.sample]["runs"]:
        run_fastqs.append(f"{config['fastq_dir']}/runs/{run}{wildcards.suffix}")

    return run_fastqs


public_samples = [sample for sample, values in SAMPLEDICT.items() if "runs" in values]
keep_fastqs = (WORKFLOW == "download_fastq") or config.get("keep_downloaded_fastq")


rule runs2sample:
    """
    Concatenate a single run or multiple runs together into a fastq
    """
    input:
        get_runs_from_sample,
    output:
             expand("{fastq_dir}/{{sample}}{{suffix}}", **config) if keep_fastqs else \
        temp(expand("{fastq_dir}/{{sample}}{{suffix}}", **config))
    log:
        expand("{log_dir}/run2sample/{{sample}}{{suffix}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/run2sample/{{sample}}{{suffix}}.benchmark.txt", **config)[0]
    wildcard_constraints:
        sample="|".join(public_samples) if len(public_samples) > 0 else "$a",
    run:
        shell("cp {input[0]} {output}")

        # now append all the later ones
        for i in range(1, len(input)):
            inputfile = input[i]
            shell("cat {inputfile} >> {output}")
