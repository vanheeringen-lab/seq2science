rule fastqc:
    input:
        f"{{path}}/{{fname}}.{config['fqext']}.gz"
    output:
        "{path}/{fname}_fastqc.html",
        "{path}/{fname}_fastqc.zip"
    conda:
        "../envs/qc.yaml"
    shell:
        "fastqc {input} -O {wildcards.path} --quiet"
