rule id2sra:
    # easy ascp install:
    # https://gist.github.com/mfansler/71f09c8b6c9a95ec4e759a8ffc488be3
    output:
        temp(directory(expand("{fastq_dir}/sra/{{sample}}", **config)))
    log:
        expand("logs/id2sra/{{sample}}.log", **config)
    threads: 1
    resources:
        parallel_downloads=1
    wildcard_constraints:
        sample="(GSM|SRR|ERR|DRR)\d+"
    conda:
        "../envs/get_fastq.yaml"
    shell:
        """
        if [[ {wildcards.sample} =~ GSM ]]; then
            IDS=$(esearch -db sra -query {wildcards.sample} | efetch --format runinfo | cut -d ',' -f 1 | grep SRR);
            TYPE_U=SRR;
        else
            IDS={wildcards.sample};
            TYPE_U=$(echo {wildcards.sample} | grep 'SRR|ERR|DRR' -E -o);
        fi;

        ASCPPATH=$(which ascp); KEYPATH=$(find / -name asperaweb_id_dsa.openssh -mount 2>/dev/null) || true;
        TYPE_L=$(echo $TYPE_U | tr '[:upper:]' '[:lower:]');
        for ID in $IDS;
        do
            PREFIX=$(echo $ID | cut -c1-6);
            SUFFIX=$(echo $ID | cut -c10-99 | xargs printf "%03d\n");

            URL_EU=\"era-fasp@fasp.sra.ebi.ac.uk:/vol1/$TYPE_L/$PREFIX/$SUFFIX/$ID\";
            URL_US=\"anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$TYPE_U/$PREFIX/$ID/$ID.sra\";
            WGET_URL=$(esearch -db sra -query {wildcards.sample} | efetch --format runinfo | cut -d ',' -f 10 | grep http);
            # first try European server, if not successful, try NCBI
            $ASCPPATH -i $KEYPATH -P33001 -T -d -k 0 -Q -l 2G -m 250M $URL_EU {output[0]}  > {log} 2>&1 ||
            $ASCPPATH -i $KEYPATH         -T -d -k 0 -Q -l 2G -m 250M $URL_US {output[0]} >> {log} 2>&1 || 
            # if both ascp servers do not work, then simple wget
            wget -O {output[0]}/{wildcards.sample} -a {log} -nv $WGET_URL 
        done;
        """


rule sra2fastq:
    # https://www.biostars.org/p/111040/
    input:
        expand("{fastq_dir}/sra/{{sample}}", fastq_dir=config['fastq_dir'])
    output:
        expand("{fastq_dir}/fastq/{{sample}}.fastq.gz", **config)
    log:
        expand("logs/sra2fastq/{{sample}}.log", **config)
    threads: 8
    params:
        files=expand("{fastq_dir}/sra/{{sample}}/*", **config)[0]
    conda:
        "../envs/get_fastq.yaml"
    shell:
        f"""
        parallel-fastq-dump -s {{params.files}} -O {config['fastq_dir']}/fastq --split-spot --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --threads {{threads}} --gzip >> {{log}} 2>&1;
        for base in $(basename -s .sra -a {{params.files}}); do
            file=$"{config['fastq_dir']}/fastq/"$base".fastq.gz"
            if [[ "$file" != {{output[0]}} ]]; then
                cat $file >> {{output[0]}} && rm $file;
            fi;
        done
        """
