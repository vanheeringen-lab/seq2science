# TODO: ASCPPATH & KEYPATH
rule gsm2sra:
    output:
        temp(directory(expand("{result_dir}/sra/{{sample}}", **config)))
    log:
        expand("{log_dir}/gsm2sra/{{sample}}.log", **config)
    resources:
        parallel_downloads=1
    wildcard_constraints:
        sample="(GSM|SRR|DRR)\d+"
    conda:
        "../envs/get_fastq.yaml"
    shell:
        """
        ASCPPATH=$HOME/.aspera/connect/bin/ascp;
        KEYPATH=$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh;

        if [[ {wildcards.sample} =~ GSM ]]; then
            IDS=$(esearch -db sra -query {wildcards.sample} | efetch --format runinfo | cut -d ',' -f 1 | grep SRR);
            TYPE_U=SRR;
        else
            IDS={wildcards.sample};
            TYPE_U=$(echo {wildcards.sample} | grep 'SRR|ERR|DRR' -E -o);
        fi;
        TYPE_L=$(echo $TYPE_U | tr '[:upper:]' '[:lower:]');

        for ID in $IDS;
        do
            PREFIX=$(echo $ID | cut -c1-6);
            SUFFIX=$(echo $ID | cut -c10-99 | xargs printf "%03d\n");

            URL_ENA1="era-fasp@fasp.sra.ebi.ac.uk:/vol1/$TYPE_L/$PREFIX/$SUFFIX/$ID";
            URL_ENA2="era-fasp@fasp.sra.ebi.ac.uk:/vol1/$TYPE_L/$PREFIX/$ID";
            URL_NCBI="anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$TYPE_U/$PREFIX/$ID/$ID.sra";
            WGET_URL=$(esearch -db sra -query {wildcards.sample} | efetch --format runinfo | cut -d ',' -f 10 | grep http);

            # first try the ENA (which has at least two different filing systems), if not successful, try NCBI
            $ASCPPATH -i $KEYPATH -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_ENA1 {output[0]}  > {log} 2>&1 ||
            $ASCPPATH -i $KEYPATH -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_ENA2 {output[0]}  > {log} 2>&1 ||
            $ASCPPATH -i $KEYPATH         -T -d -k 0 -Q -l 2G -m 250M $URL_NCBI {output[0]}  > {log} 2>&1 ||
            # if none of the ascp servers work, then simply wget
            wget -O {output[0]}/{wildcards.sample} -a {log} -nv $WGET_URL 
        done;
        """


rule sra2fastqPE_split:
    input:
        rules.gsm2sra.output
    output:
        expand("{fastq_dir}/fastq/{{sample}}_{fqext1}.{fqsuffix}.gz", **config),
        expand("{fastq_dir}/fastq/{{sample}}_{fqext2}.{fqsuffix}.gz", **config)
    log:
        expand("{log_dir}/sra2fastq/{{sample}}.log", **config)
    threads: 8
    params:
        files=expand("{result_dir}/sra/{{sample}}/*", **config)[0],
        params=config['PE_split']
    conda:
        "../envs/get_fastq.yaml"
    shell:
        f"""
        parallel-fastq-dump -s {{params.files}} -O {config['result_dir']}/fastq {{params.params}} --threads {{threads}} --gzip >> {{log}} 2>&1

        # rename the SRRs to GSMs
        GSM=$(basename {{input}})
        SRR=$(basename {{params.files}})
        FILES={config['result_dir']}/*.{config['fqsuffix']}.gz
        rcmd=s/$SRR/$GSM/
        rename -v $rcmd $FILES >> {{log}} 2>&1
        """


rule sra2fastqSE_splot:
    input:
        rules.gsm2sra.output
    output:
        expand("{result_dir}/fastq/{{sample}}.{fqsuffix}.gz", **config)
    log:
        expand("{log_dir}/sra2fastq/{{sample}}.log", **config)
    threads: 40
    params:
        files=expand("{result_dir}/sra/{{sample}}/*", **config)[0],
        params=config['SE_splot']
    conda:
        "../envs/get_fastq.yaml"
    shell:
        f"""
        parallel-fastq-dump -s {{params.files}} -O {config['result_dir']}/fastq {{params.params}} --threads {{threads}} --gzip >> {{log}} 2>&1
        
        # rename the SRR to GSM
        GSM=$(basename {{input}})
        SRR=$(basename {{params.files}})
        src='{config['result_dir']}/fastq/'$SRR'_*.{config['fqsuffix']}.gz'
        dst='{config['result_dir']}/fastq/'$GSM'.{config['fqsuffix']}.gz'
        echo $src; echo $dst;
        mv -v $src $dst >> {{log}} 2>&1
        """
