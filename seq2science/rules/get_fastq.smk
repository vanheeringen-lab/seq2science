import glob
import os
import re
import time


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
    wildcard_constraints:
        sample="(GSM|SRR|ERR|DRR)\d+",
    params:
        ascp_path=config.get("ascp_path", "NO_ASCP_PATH_PROVIDED"),
        ascp_key=config.get("ascp_key", "NO_ASCP_key_PROVIDED"),
    run:
        with FileLock(layout_cachefile_lock):
            shell(
                """
        echo "starting lookup of the sample in the sra database:" > {log}
        if [[ {wildcards.sample} =~ GSM ]]; then
            IDS=$(esearch -db sra -query {wildcards.sample} | efetch --format runinfo | cut -d ',' -f 1 | grep SRR);
            TYPE_U=SRR;
        else
            IDS={wildcards.sample};
            TYPE_U=$(echo {wildcards.sample} | grep 'SRR|ERR|DRR' -E -o);
        fi;
        TYPE_L=$(echo $TYPE_U | tr '[:upper:]' '[:lower:]');
        echo "ids: $IDS, type (uppercase): $TYPE_U, type (lowercase): $TYPE_L" >> {log}

        declare -a URLS
        for ID in $IDS;
        do
            sleep 2s;
            WGET_URL=$(esearch -db sra -query $ID | efetch --format runinfo | grep $ID | cut -d ',' -f 10 | grep http);
            URLS+=($WGET_URL);
        done;

        mkdir {output}
        printf '%s\n' "${{IDS[@]}}" > {output}/ids
        printf '%s\n' "${{URLS[@]}}" > {output}/urls
        echo $TYPE_L > {output}/type_l
        echo $TYPE_U > {output}/type_u
        """
            )
            time.sleep(2)

        shell(
            """
        read TYPE_L < {output}/type_l
        read TYPE_U < {output}/type_u
        readarray -t IDS < {output}/ids
        readarray -t URLS < {output}/urls
        rm {output}/*
        for i in ${{!IDS[@]}}; do
            ID="${{IDS[i]}}"

            PREFIX=$(echo $ID | cut -c1-6);
            SUFFIX=$(echo -n $ID | tail -c 3);

            URL_ENA1="era-fasp@fasp.sra.ebi.ac.uk:/vol1/$TYPE_L/$PREFIX/$SUFFIX/$ID";
            URL_ENA2="era-fasp@fasp.sra.ebi.ac.uk:/vol1/$TYPE_L/$PREFIX/$ID";
            URL_NCBI="anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$TYPE_U/$PREFIX/$ID/$ID.sra";
            WGET_URL=${{URLS[i]}}

            echo "trying to download $ID from, respectively: \n$URL_ENA1 \n$URL_ENA2 \n$URL_NCBI \n$WGET_URL" >> {log}
            # first try the ENA (which has at least two different filing systems), if not successful, try NCBI
            # if none of the ascp servers work, or ascp is not defined in the config, then simply wget
            {params.ascp_path} -i {params.ascp_key} -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_ENA1 {output[0]} >> {log} 2>&1 ||
            {params.ascp_path} -i {params.ascp_key} -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_ENA2 {output[0]} >> {log} 2>&1 ||
            {params.ascp_path} -i {params.ascp_key}         -T -d -k 0 -Q -l 2G -m 250M $URL_NCBI {output[0]} >> {log} 2>&1 ||
            (mkdir -p {output[0]} >> {log} 2>&1 && wget -O {output[0]}/$ID -a {log} -nv $WGET_URL >> {log} 2>&1)
        done;
        """
        )



rule sra2fastq_SE:
    """
    Downloaded (raw) SRAs are converted to single-end fastq files.
    """
    input:
        rules.id2sra.output,
    output:
        expand("{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config),
    log:
        expand("{log_dir}/sra2fastq_SE/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sra2fastq_SE/{{sample}}.benchmark.txt", **config)[0]
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        """
        # setup tmp dir
        tmpdir={config[sra_dir]}/tmp/{wildcards.sample}
        mkdir -p $tmpdir; trap "rm -rf $tmpdir" EXIT

        # dump to tmp dir
        parallel-fastq-dump -s {input}/* -O $tmpdir {config[splot]} \
        --threads {threads} --gzip >> {log} 2>&1

        # rename file and move to output dir
        for f in $(ls -1q $tmpdir | grep -oP "^[^_]+" | uniq); do
            dst={config[fastq_dir]}/{wildcards.sample}.{config[fqsuffix]}.gz
            cat "${{tmpdir}}/${{f}}_pass.fastq.gz" >> $dst
        done
        """


ruleorder: renamefastq_PE > sra2fastq_PE


rule sra2fastq_PE:
    """
    Downloaded (raw) SRAs are converted to paired-end fastq files.
    Forward and reverse samples will be switched if forward/reverse names are not lexicographically ordered.
    """
    input:
        rules.id2sra.output,
    output:
        expand("{fastq_dir}/{{sample}}_{fqext}.{fqsuffix}.gz", **config),
    log:
        expand("{log_dir}/sra2fastq_PE/{{sample}}.log", **config),
    benchmark:
        expand("{benchmark_dir}/sra2fastq_PE/{{sample}}.benchmark.txt", **config)[0]
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        """
        # setup tmp dir
        tmpdir={config[sra_dir]}/tmp/{wildcards.sample}
        mkdir -p $tmpdir; trap "rm -rf $tmpdir" EXIT

        # dump to tmp dir
        parallel-fastq-dump -s {input}/* -O $tmpdir {config[split]} \
        --threads {threads} --gzip >> {log} 2>&1

        # rename files and move to output dir
        for f in $(ls -1q $tmpdir | grep -oP "^[^_]+" | uniq); do
            dst_1={config[fastq_dir]}/{wildcards.sample}_{config[fqext1]}.{config[fqsuffix]}.gz
            dst_2={config[fastq_dir]}/{wildcards.sample}_{config[fqext2]}.{config[fqsuffix]}.gz
            cat "${{tmpdir}}/${{f}}_pass_1.fastq.gz" >> $dst_1
            cat "${{tmpdir}}/${{f}}_pass_2.fastq.gz" >> $dst_2
        done
        """


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
