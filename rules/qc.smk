rule fastqc:
    input:
        "{path}/{fname}.fastq.gz"
    output:
        "{path}/{fname}_fastqc.html",
        "{path}/{fname}_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    shell:
        "fastqc {input} -O {wildcards.path} --quiet"
