

rule gimmemotifs:
    input:
        pass
    output:
        pass
    conda:
        "../envs/gimme.yaml"
    shell:
        """
        
        """


rule plot_gimmemotifs:
    input:
        rules.gimmemotifs.output
    output:
        expand("{qc_dir}/gimmemotifs/{{assembly}}-{{peak_caller}}_mqc.png", **config)
    conda:
        "../envs/gimme.yaml"
    script:
        pass