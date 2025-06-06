"""
all rules/logic related downloading public fastqs should be here.
"""
localrules: run2sra, ena2fastq_SE, ena2fastq_PE, gsa_or_encode2fastq_SE, gsa_or_encode2fastq_PE, runs2sample


import os


rule run2sra:
    """
    Download the SRA of a sample by its unique identifier.
    """
    output:
        temp(expand("{sra_dir}/{{run}}/{{run}}/{{run}}.sra", **config)),
    log:
        expand("{log_dir}/run2sra/{{run}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/run2sra/{{run}}.benchmark.txt", **config)[0]
    message: EXPLAIN["run2sra"]
    resources:
        parallel_downloads=1,
    params:
        outdir=lambda wildcards: f"{config['sra_dir']}/{wildcards.run}",
    conda:
        "../envs/get_fastq.yaml"
    wildcard_constraints:
        run="[DES]RR\d+",
    retries: 2
    shell:
        """
        # move to output dir since somehow prefetch sometimes puts files in the cwd...
        # and remove the top level folder since prefetch will assume we are done otherwise
        mkdir -p {params.outdir}; cd {params.outdir}; rm -r {wildcards.run}

        # TODO: for loop can be removed if we dont see the debug message
        # three attempts
        for i in {{1..3}}
        do
            # acquire a lock
            (
                flock -w 30 200 || continue
                sleep 2
            ) 200>{PYSRADB_CACHE_LOCK}

            # dump
            prefetch --max-size u --output-directory ./ --log-level debug --progress {wildcards.run} \
            >> {log} 2>&1 && break
            
            # retry
            echo "DEBUG: prefetch try ${{i}} of 3 failed" >> {log}
            sleep 10
        done

        # TODO: section can be removed if we dont see the debug message
        # bug report: https://github.com/ncbi/sra-tools/issues/533
        if [[ -f "{params.outdir}/{wildcards.run}.sra" ]]; then
            echo "DEBUG: moving output to correct directory" >> {log}
            mkdir -p {params.outdir}/{wildcards.run}
            mv {params.outdir}/{wildcards.run}.sra {output}
        fi

        # TODO: section can be removed if we dont see the debug message
        # If an sralite file was downloaded instead of a sra file, just rename it
        if [[ -f "{params.outdir}/{wildcards.run}.sralite" ]]; then
            echo "DEBUG: renaming SRAlite" >> {log}
            mkdir -p {params.outdir}/{wildcards.run}
            mv {params.outdir}/{wildcards.run}.sralite {output}
        fi
        """


# ENA > SRA
ruleorder: ena2fastq_SE > gsa_or_encode2fastq_SE > sra2fastq_SE
ruleorder: ena2fastq_PE > gsa_or_encode2fastq_PE > sra2fastq_PE
# PE > SE
ruleorder: ena2fastq_PE > ena2fastq_SE
ruleorder: sra2fastq_PE > sra2fastq_SE
ruleorder:gsa_or_encode2fastq_PE > gsa_or_encode2fastq_SE


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
        run="|".join(SRA_SINGLE_END) if len(SRA_SINGLE_END) else "$a",  # only try to dump (single-end) SRA samples
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    priority: 1
    retries: 2
    shell:
        """
        # move to output dir since somehow parallel-fastq-dump sometimes puts files in the cwd...
        mkdir -p {output.tmpdir}; cd {output.tmpdir}

        # acquire the lock
        (
            flock -w 30 200 || exit 1
            sleep 3
        ) 200>{PYSRADB_CACHE_LOCK}

        # dump to tmp dir
        parallel-fastq-dump -s {input} -O {output.tmpdir} \
        --threads {threads} --split-spot --skip-technical --dumpbase --readids \
        --clip --read-filter pass --defline-seq '@$ac.$si.$sg/$ri' \
        --defline-qual '+' --gzip >> {log} 2>&1

        # rename file and move to output dir
        mv {output.tmpdir}/*.fastq.gz {output.fastq}
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
        run="|".join(SRA_PAIRED_END) if len(SRA_PAIRED_END) else "$a",  # only try to dump (paired-end) SRA samples
    conda:
        "../envs/get_fastq.yaml"
    priority: 1
    retries: 2
    shell:
        """
        # move to output dir since somehow parallel-fastq-dump sometimes puts files in the cwd...
        mkdir -p {output.tmpdir}; cd {output.tmpdir}

        # acquire the lock
        (
            flock -w 30 200 || exit 1
            sleep 3
        ) 200>{PYSRADB_CACHE_LOCK}

        # dump to tmp dir
        parallel-fastq-dump -s {input} -O {output.tmpdir} \
        --threads {threads} --split-3 --skip-technical --dumpbase \
        --readids --clip --read-filter pass --defline-seq '@$ac.$si.$sg/$ri' \
        --defline-qual '+' --gzip >> {log} 2>&1
        
        # check if the files exist
        if ! compgen -G "{output.tmpdir}/*_1*" > /dev/null ; then printf "ERROR: Couldn't find read 1.fastq after dumping! Perhaps this is not a paired-end file?\n" >> {log} 2>&1; fi
        if ! compgen -G "{output.tmpdir}/*_2*" > /dev/null ; then printf "ERROR: Couldn't find read 2.fastq after dumping! Perhaps this is not a paired-end file?\n" >> {log} 2>&1; fi

        # rename file and move to output dir
        mv {output.tmpdir}/*_1* {output.fastq[0]} >> {log} 2>&1
        mv {output.tmpdir}/*_2* {output.fastq[1]} >> {log} 2>&1
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
    priority: 1
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
    priority: 1
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


rule gsa_or_encode2fastq_SE:
    """
    Download single-end fastq files from the GSA or ENCODE DCC.
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
        run="|".join(GSA_OR_ENCODE_SINGLE_END) if len(GSA_OR_ENCODE_SINGLE_END) else "$a",
    priority: 1
    retries: 2
    run:
        shell("mkdir -p {config[fastq_dir]}/tmp/ >> {log} 2>&1")
        url = RUN2DOWNLOAD[wildcards.run]

        shell("wget {url} -O {output} --waitretry 20 >> {log} 2>&1")

        if not (os.path.exists(output[0]) and os.path.getsize(output[0]) > 0):
            shell('echo "Something went wrong.. The downloaded file was empty!" >> {log} 2>&1')
            shell("exit 1 >> {log} 2>&1")


rule gsa_or_encode2fastq_PE:
    """
    Download paired-end fastq files from the GSA or ENCODE DCC.
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
        run="|".join(GSA_OR_ENCODE_PAIRED_END) if len(GSA_OR_ENCODE_PAIRED_END) else "$a",
    priority: 1
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
    priority: 2
    wildcard_constraints:
        sample="|".join(public_samples) if len(public_samples) > 0 else "$a",
    run:
        shell("cp {input[0]} {output}")

        # now append all the later ones
        for i in range(1, len(input)):
            inputfile = input[i]
            shell("cat {inputfile} >> {output}")
