rule id2sra:
    output:
        temp(directory(expand("{result_dir}/{sra_dir}/{{sample}}", **config)))
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
        ascp_path=config['ascp_path'],
        ascp_key= config['ascp_key']
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
            WGET_URL=$(esearch -db sra -query {wildcards.sample} | efetch --format runinfo | cut -d ',' -f 10 | grep http);

            echo "trying to download $ID from, respectively: \n$URL_ENA1 \n$URL_ENA2 \n$URL_NCBI \n$WGET_URL" >> {log}
            # first try the ENA (which has at least two different filing systems), if not successful, try NCBI
            # if none of the ascp servers work, or ascp is not defined in the config, then simply wget
            {params.ascp_path} -i {params.ascp_key} -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_ENA1 {output[0]} >> {log} 2>&1 ||
            {params.ascp_path} -i {params.ascp_key} -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_ENA2 {output[0]} >> {log} 2>&1 ||
            {params.ascp_path} -i {params.ascp_key}         -T -d -k 0 -Q -l 2G -m 250M $URL_NCBI {output[0]} >> {log} 2>&1 ||
            (mkdir {output[0]} >> {log} 2>&1 && wget -O {output[0]}/{wildcards.sample} -a {log} -nv $WGET_URL >> {log} 2>&1)
        done;
        """


rule sra2fastq_SE:
    input:
        rules.id2sra.output
    output:
        expand("{result_dir}/{fastq_dir}/{{sample}}.{fqsuffix}.gz", **config)
    log:
        expand("{log_dir}/sra2fastq_SE/{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/sra2fastq_SE/{{sample}}.benchmark.txt", **config)[0]
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        f"""
        parallel-fastq-dump -s {{input}}/* -O {config['result_dir']}/{config['fastq_dir']} {config['splot']} --threads {{threads}} --gzip >> {{log}} 2>&1

        # rename the SRR to GSM
        GSM=$(basename {{input}})
        SRR=$(basename {{input}}/*)
        src='{config['result_dir']}/{config['fastq_dir']}/'$SRR'_pass.{config['fqsuffix']}.gz'
        dst='{config['result_dir']}/{config['fastq_dir']}/'$GSM'.{config['fqsuffix']}.gz'
        mv -v $src $dst >> {{log}} 2>&1
        """


rule sra2fastq_PE:
    input:
        rules.id2sra.output
    output:
        expand("{result_dir}/{fastq_dir}/{{sample}}_{fqext}.{fqsuffix}.gz", **config)
    log:
        expand("{log_dir}/sra2fastq_PE/{{sample}}.log", **config)
    benchmark:
        expand("{benchmark_dir}/sra2fastq_PE/{{sample}}.benchmark.txt", **config)[0]
    threads: 8
    conda:
        "../envs/get_fastq.yaml"
    shell:
        f"""
        parallel-fastq-dump -s {{input}}/* -O {config['result_dir']}/{config['fastq_dir']} {config['split']} --threads {{threads}} --gzip >> {{log}} 2>&1

        # rename the SRRs to GSMs
        SRR=$(basename {{input}}/*)
        GSM=$(basename {{input}})
        file1="$(find "$(dirname {{output[0]}})" -name "$SRR*.fastq.gz" | sort | head -n 1)"
        mv "$file1" "$(echo "$file1" | sed -e "s/$SRR/$GSM/" -e 's/.fastq/.{config['fqsuffix']}/' -e 's/pass_1/{config['fqext1']}/')"
        file1="$(find "$(dirname {{output[0]}})" -name "$SRR*.fastq.gz" | sort | tail -n 1)"
        mv "$file1" "$(echo "$file1" | sed -e "s/$SRR/$GSM/" -e 's/.fastq/.{config['fqsuffix']}/' -e 's/pass_2/{config['fqext2']}/')"
        """
