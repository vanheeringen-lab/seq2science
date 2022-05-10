"""
all rules/logic related to making gene count tables should be here.
"""

def get_counts(wildcards):
    """return a list with salmon transcript abundances for all samples in an assembly"""
    quant_dirs = []
    for sample in treps[treps["assembly"] == ORI_ASSEMBLIES[wildcards.assembly]].index:
        quant_dirs.append(f"{{result_dir}}/salmon/{wildcards.assembly}-{sample}/quant.sf")
    return expand(quant_dirs, **config)

def get_names(wildcards):
    """return the descriptive>technical_replicate>sample names of the samples"""
    names = []
    for sample in treps[treps["assembly"] == ORI_ASSEMBLIES[wildcards.assembly]].index:
        names.append(rep_to_descriptive(sample))
    return names


if config["quantifier"] == "salmon" and config["tpm2counts"] == "tximeta":

    rule linked_txome:
        """
        Generate a linked transcriptome for tximeta

        Also creates a symlink to the gtf in an Ensembl format (required by tximeta)

        Required to converting salmon output (estimated transcript abundances) to gene counts
        """
        input:
            fasta=rules.get_transcripts.output,
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf", **config),
            index_dir=rules.salmon_index.output,
        output:
            index=expand(f"{{genome_dir}}/{{{{assembly}}}}/index/tximeta/{salmon_index_output}_txome.json", **config),
        params:
            source=config["txi_source"],
            organism=config["txi_organism"],
            release=config["txi_release"],
        log:
            expand("{log_dir}/quantification/{{assembly}}-linked_txome.log", **config),
        conda: "../envs/tximeta.yaml"
        resources: R_scripts=1  # conda's R can have issues when starting multiple times
        script:
            f"{config['rule_dir']}/../scripts/generate_linked_txome.R"

    rule txi_count_matrix:
        """
        Convert transcript abundance estimations to gene count estimations and merge
        gene counts per assembly.

        Also outputs a single cell experiment object similar to ARMOR
        (https://github.com/csoneson/ARMOR).

        Only works with Ensembl assemblies.
        """
        input:
            cts=get_counts,
            linked_txome=rules.linked_txome.output,
        output:
            counts=expand("{counts_dir}/{{assembly}}-counts.tsv", **config),
            tpms=expand("{counts_dir}/{{assembly}}-TPM.tsv",**config),
            lengths=expand("{counts_dir}/{{assembly}}-gene_lengths.tsv",**config),
            SCE=expand("{counts_dir}/{{assembly}}-se.rds", **config),
        log:
            expand("{log_dir}/quantification/{{assembly}}-txi_counts_matrix.log", **config),
        params:
            reps=lambda wildcards, input: input,  # help resolve changes in input files
            names=lambda wildcards: get_names(wildcards),
        message: EXPLAIN["count_matrix_txi"]
        conda: "../envs/tximeta.yaml"
        resources: R_scripts=1  # conda's R can have issues when starting multiple times
        script:
            f"{config['rule_dir']}/../scripts/quant_to_counts.R"


elif config["quantifier"] == "salmon" and config["tpm2counts"] == "pytxi":

    rule pytxi_count_matrix:
        """
        Convert transcript abundance estimations to gene count estimations and merge 
        gene counts per assembly.

        If the GTF is not used, a taxonomy ID is required for myGene.info.
        Seq2science will try to read this ID from the genomepy README.txt or assembly database. 
        """
        input:
            cts=get_counts,
            fa=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf",**config),
        output:
            counts=expand("{counts_dir}/{{assembly}}-counts.tsv",**config),
            tpms=expand("{counts_dir}/{{assembly}}-TPM.tsv",**config),
            lengths=expand("{counts_dir}/{{assembly}}-gene_lengths.tsv",**config),
        params:
            reps=lambda wildcards, input: input,  # help resolve changes in input files
            names=lambda wildcards: get_names(wildcards),
            from_gtf=config["tx2gene_from_gtf"],
        log:
            expand("{log_dir}/quantification/{{assembly}}-pytxi_counts_matrix.log",**config),
        message: EXPLAIN["count_matrix_pytxi"]
        conda: "../envs/pytxi.yaml"
        script:
            f"{config['rule_dir']}/../scripts/pytxi.py"

else:

    def get_counts(wildcards):
        count_tables = []
        for sample in treps[treps["assembly"] == ORI_ASSEMBLIES[wildcards.assembly]].index:
            count_tables.append(f"{{counts_dir}}/{wildcards.assembly}-{sample}.counts.tsv")
        return expand(count_tables, **config)

    rule count_matrix:
        """
        Combine count tables into one count matrix per assembly
        """
        input:
            cts=get_counts,
        output:
            expand("{counts_dir}/{{assembly}}-counts.tsv", **config),
        params:
            reps=lambda wildcards, input: input,  # help resolve changes in input files
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-counts_matrix.log", **config),
        run:
            import pandas as pd
            import sys

            sys.stdout = open(log[0],'w')
            sys.stderr = sys.stdout

            counts = pd.DataFrame()
            for sample in input.cts:
                col = pd.read_csv(
                    sample,
                    sep="\t",
                    index_col=0,
                    header=None,
                    skiprows=2 if config["quantifier"] == "featurecounts" else 0,
                    usecols=[0, 6] if config["quantifier"] == "featurecounts" else [0, 1],
                    skipfooter=5 if config["quantifier"] == "htseq" else 0,
                )
                sample_name = sample.split(wildcards.assembly + "-")[1].split(".counts.tsv")[0]
                col.columns = [rep_to_descriptive(sample_name)]
                counts = pd.concat([counts, col], axis=1)

            counts.index.name = "gene"
            counts.to_csv(output[0], sep="\t")

    rule tpm_matrix:
        """
        Create a TPM normalized table from a counts table.
        """
        input:
            cts=rules.count_matrix.output,
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.annotation.gtf",**config),
        output:
            tpms = expand("{counts_dir}/{{assembly}}-TPM.tsv",**config),
            lengths = expand("{counts_dir}/{{assembly}}-gene_lengths.tsv",**config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-tpm_matrix.log",**config),
        message: EXPLAIN["tpm_matrix"]
        script:
            f"{config['rule_dir']}/../scripts/counts2tpm.py"


if config.get("dexseq"):

    def get_dexseq_counts(wildcards):
        count_tables = []
        for sample in treps[treps["assembly"] == ORI_ASSEMBLIES[wildcards.assembly]].index:
            count_tables.append(f"{{counts_dir}}/{wildcards.assembly}-{sample}.DEXSeq_counts.tsv")
        return expand(count_tables, **config)

    rule dexseq_count_matrix:
        """
        Combine DEXSeq counts into one count matrix per assembly
        for use in function `DEXSeqDataSet()`
        """
        input:
            cts=get_dexseq_counts,
        output:
            expand("{counts_dir}/{{assembly}}-DEXSeq_counts.tsv", **config),
        params:
            reps=lambda wildcards, input: input,  # help resolve changes in input files
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-DEXSeq_counts_matrix.log", **config),
        run:
            import pandas as pd
            import sys

            sys.stdout = open(log[0],'w')
            sys.stderr = sys.stdout

            counts = pd.DataFrame()
            for sample in input.cts:
                col = pd.read_csv(sample, sep="\t", index_col=0, header=None, skipfooter=5)
                sample_name = sample.split(wildcards.assembly + "-")[1].split(".DEXSeq_counts.tsv")[0]
                col.columns = [rep_to_descriptive(sample_name)]
                counts = pd.concat([counts, col], axis=1)

            counts.index.name = "exon"
            counts.to_csv(output[0], sep="\t")
