rule id2sra:
    """
    Download the SRA of a sample by its unique identifier.

    Tries first downloading with the faster ascp protocol, if that fails it falls back on the slower http protocol.
    """
    output:
        temp(directory(expand("{sra_dir}/{{sample}}", **config)))
    log:
        expand("{log_dir}/id2sra/{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/id2sra/{{sample}}.benchmark.txt", **config)[0]
    resources:
        parallel_downloads=1
    wildcard_constraints:
        sample="(GSM|SRR|ERR|DRR)\d+"
    conda:
        "../envs/get_fastq.yaml"
    params:
        ascp_path=config.get('ascp_path', "NO_ASCP_PATH_PROVIDED"),
        ascp_key= config.get('ascp_key', "NO_ASCP_key_PROVIDED")
    shell:
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

        for ID in $IDS;
        do
            PREFIX=$(echo $ID | cut -c1-6);
            SUFFIX=$(echo $ID | cut -c10-99 | xargs printf "%03d\n");

            URL_ENA1="era-fasp@fasp.sra.ebi.ac.uk:/vol1/$TYPE_L/$PREFIX/$SUFFIX/$ID";
            URL_ENA2="era-fasp@fasp.sra.ebi.ac.uk:/vol1/$TYPE_L/$PREFIX/$ID";
            URL_NCBI="anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$TYPE_U/$PREFIX/$ID/$ID.sra";
            WGET_URL=$(esearch -db sra -query $ID | efetch --format runinfo | grep $ID | cut -d ',' -f 10 | grep http);

            echo "trying to download $ID from, respectively: \n$URL_ENA1 \n$URL_ENA2 \n$URL_NCBI \n$WGET_URL" >> {log}
            # first try the ENA (which has at least two different filing systems), if not successful, try NCBI
            # if none of the ascp servers work, or ascp is not defined in the config, then simply wget
            {params.ascp_path} -i {params.ascp_key} -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_ENA1 {output[0]} >> {log} 2>&1 ||
            {params.ascp_path} -i {params.ascp_key} -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_ENA2 {output[0]} >> {log} 2>&1 ||
            {params.ascp_path} -i {params.ascp_key}         -T -d -k 0 -Q -l 2G -m 250M $URL_NCBI {output[0]} >> {log} 2>&1 ||
            (mkdir -p {output[0]} >> {log} 2>&1 && wget -O {output[0]}/$ID -a {log} -nv $WGET_URL >> {log} 2>&1)
        done;

        #if the folder contains multiple files, concatenate them
        if [[ $(ls -1q {output[0]} | wc -l) > 1 ]]; then
          catfile={output[0]}/$TYPE_U'_concatenated'
          for file in {output[0]}/*; do
            cat $file >> $catfile && rm $file
          done
        fi
        """


rule sra2fastq_SE:
    """
    Downloaded (raw) SRAs are converted to single-end fastq files.
    """
    input:
        rules.id2sra.output
    output:
        expand("{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config)
    log:
        expand("{log_dir}/sra2fastq_SE/{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/sra2fastq_SE/{{sample}}.benchmark.txt", **config)[0]
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        f"""
        parallel-fastq-dump -s {{input}}/* -O {config['fastq_dir']} {config['splot']} --threads {{threads}} --gzip >> {{log}} 2>&1

        # rename the SRR to GSM
        SRR=$(basename {{input}}/*)
        GSM=$(basename {{input}})
        src='{config['fastq_dir']}/'$SRR'_pass.{config['fqsuffix']}.gz'
        dst='{config['fastq_dir']}/'$GSM'.{config['fqsuffix']}.gz'
        if [[ ! $src = $dst ]]; then
              mv -v $src $dst >> {{log}} 2>&1
        fi
        """


ruleorder: renamefastq_PE > sra2fastq_PE
rule sra2fastq_PE:
    """
    Downloaded (raw) SRAs are converted to paired-end fastq files.
    Forward and reverse samples will be switched if forward/reverse names are not lexicographically ordered.
    """
    input:
        rules.id2sra.output
    output:
        expand("{fastq_dir}/{{sample}}_{fqext}.{fqsuffix}.gz", **config)
    log:
        expand("{log_dir}/sra2fastq_PE/{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/sra2fastq_PE/{{sample}}.benchmark.txt", **config)[0]
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        f"""
        parallel-fastq-dump -s {{input}}/* -O {config['fastq_dir']} {config['split']} --threads {{threads}} --gzip >> {{log}} 2>&1

        # renaming SRRs to GSMs
        SRR=$(basename {{input}}/*)
        GSM=$(basename {{input}})
        if [[ ! $SRR = $GSM ]]; then
          echo "\nrenaming SRRs to GSMs:" >> {{log}}

          src='{config['fastq_dir']}/'$SRR'_{config['fqext1']}.{config['fqsuffix']}.gz'
          dst='{config['fastq_dir']}/'$GSM'_{config['fqext1']}.{config['fqsuffix']}.gz'
          mv -v $src $dst >> {{log}} 2>&1

          src='{config['fastq_dir']}/'$SRR'_{config['fqext2']}.{config['fqsuffix']}.gz'
          dst='{config['fastq_dir']}/'$GSM'_{config['fqext2']}.{config['fqsuffix']}.gz'
          mv -v $src $dst >> {{log}} 2>&1
        fi

        # renaming fqexts
        all_files=$(ls -1q {config['fastq_dir']} | grep -c $GSM'_.*.gz')
        # Snakemake throws an error if variable assignment is based on a grep command without output (0 hits)
        correct_files=$(( $(ls -1q {config['fastq_dir']} | grep -c $GSM'_{config['fqext1']}.{config['fqsuffix']}.gz') + \
                          $(ls -1q {config['fastq_dir']} | grep -c $GSM'_{config['fqext2']}.{config['fqsuffix']}.gz') )) \
                      || correct_files=0
        if [[ ! $all_files = $correct_files ]]; then
            echo "\nrenaming fqexts:" >> {{log}}

            n=1
            for f in $(ls -1q {config['fastq_dir']} | grep $GSM'_.*.{config['fqsuffix']}.gz'); do
              src='{config['fastq_dir']}/'$f
              dst='{config['fastq_dir']}/'$GSM'_{config['fqext1']}.{config['fqsuffix']}.gz'
              if [[ $n = 2 ]]; then
                dst='{config['fastq_dir']}/'$GSM'_{config['fqext2']}.{config['fqsuffix']}.gz'
              fi
              mv -v $src $dst >> {{log}} 2>&1
              n=$(( $n + 1 ))
            done
        fi
        """

def get_wrong_fqext(wildcards):
    """get all samples with fqexts that are not in the config"""
    fastqs = glob.glob(os.path.join(config["fastq_dir"], wildcards.sample + "*" + config["fqsuffix"] + ".gz"))
    # exclude samples with the correct fqext
    r1 = wildcards.sample + "_" + config["fqext1"] + "." + config["fqsuffix"] + ".gz"
    r2 = wildcards.sample + "_" + config["fqext2"] + "." + config["fqsuffix"] + ".gz"
    fastqs = [f for f in fastqs if not re.match(r1+"|"+r2, f)]
    if len(fastqs) == 0:
        fastqs.append("impossible input to prevent Snakemake from running the rule without input")
    return sorted(fastqs)

rule renamefastq_PE:
    """
    Create symlinks to fastqs with incorrect fqexts (default R1/R2).
    Forward and reverse samples will be switched if forward/reverse names are not lexicographically ordered.
    """
    input:
         get_wrong_fqext
    output:
         temp(expand("{fastq_dir}/{{sample}}_{fqext}.{fqsuffix}.gz", **config))
    shell:
         """
         ln {input[0]} {output[0]}
         ln {input[1]} {output[1]}
         """