if config["quantifier"] == "salmon":

    def get_counts(wildcards):
        quant_dirs = []
        if config["quantifier"] in ["salmon", "star"]:
            for sample in treps[treps["assembly"] == wildcards.assembly].index:
                quant_dirs.append(f"{{result_dir}}/{{quantifier}}/{wildcards.assembly}-{sample}")
        return expand(quant_dirs, **config)

    def get_index(wildcards):
        index = f"{{genome_dir}}/{wildcards.assembly}/index/{{quantifier}}"
        if config["decoy_aware_index"]:
            index += "_decoy_aware"
        return expand(index, **config)

    rule linked_txome:
        """
        Generate a linked transcriptome for tximeta

        Also creates a symlink to the gtf in an Ensembl format (required by tximeta)

        Required to converting salmon output (estimated transcript abundances) to gene counts
        """
        input:
            fasta=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            index_dir=get_index,
        output:
            index=expand("{genome_dir}/{{assembly}}/index/tximeta/linked_txome.json", **config), # symlink=expand(f"{{genome_dir}}/{{{{assembly}}}}/index/tximeta/{config['tximeta']['organism']}.{{{{assembly}}}}.{config['tximeta']['release']}.gtf", **config)
        params:
            source=config["tximeta"]["source"],
            organism=config["tximeta"]["organism"],
            release=config["tximeta"]["release"],
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-linked_txome.log", **config),
        conda:
            "../envs/tximeta.yaml"
        resources:
            R_scripts=1, # conda's R can have issues when starting multiple times
        script:
            f"{config['rule_dir']}/../scripts/generate_linked_txome.R"

    rule count_matrix:
        """
        Convert transcript abundance estimations to gene count estimations and merge 
        gene counts per assembly.

        Also outputs a single cell experiment object similar to ARMOR 
        (https://github.com/csoneson/ARMOR).

        Only works with Ensembl assemblies.
        """
        input:
            linked_txome=expand("{genome_dir}/{{assembly}}/index/tximeta/linked_txome.json", **config),
            cts=get_counts,
        output:
            counts=expand("{counts_dir}/{{assembly}}-counts.tsv", **config),
            SCE=expand("{counts_dir}/{{assembly}}-se.rds", **config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-txi_counts_matrix.log", **config),
        message: explain_rule("count_matrix_txi")
        conda:
            "../envs/tximeta.yaml"
        resources:
            R_scripts=1, # conda's R can have issues when starting multiple times
        script:
            f"{config['rule_dir']}/../scripts/quant_to_counts.R"

else:

    def get_counts(wildcards):
        count_tables = []
        for sample in treps[treps["assembly"] == wildcards.assembly].index:
            count_tables.append(f"{{counts_dir}}/{wildcards.assembly}-{sample}.counts.tsv")
        return expand(count_tables, **config)

    rule count_matrix:
        """
        Combine count tables into one count matrix per assembly
        """
        input:
            cts=get_counts
        output:
            expand("{counts_dir}/{{assembly}}-counts.tsv", **config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-counts_matrix.log", **config),
        run:
            import pandas as pd
            import sys

            with open(log[0], "w") as log_file:
                sys.stderr = sys.stdout = log_file

                counts = pd.DataFrame()
                for sample in input.cts:
                    col = pd.read_csv(
                        sample,
                        sep="\t",
                        index_col=0,
                        header=None,
                        skiprows=2 if config["quantifier"] == "featurecounts" else 0,
                        usecols=[0,6] if config["quantifier"] == "featurecounts" else [0,1],
                        skipfooter=5 if config["quantifier"] == "htseq" else 0,
                    )
                    sample_name = sample.split(wildcards.assembly + "-")[1].split(".counts.tsv")[0]
                    col.columns = [sample_name]
                    counts = pd.concat([counts, col], axis=1)

                counts.index.name = "gene"
                counts.to_csv(output[0], sep="\t")
