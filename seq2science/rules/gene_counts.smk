def get_counts(wildcards):
    quant_dirs = []
    for sample in treps[treps['assembly'] == wildcards.assembly].index:
        quant_dirs.append(f"{{result_dir}}/{{quantifier}}/{wildcards.assembly}-{sample}")
    return expand(quant_dirs, **config)

if config['quantifier'] == 'salmon':
    def get_index(wildcards):
        index=f"{{genome_dir}}/{wildcards.assembly}/index/{{quantifier}}"
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
            index_dir=get_index
        output:
            index=expand("{genome_dir}/{{assembly}}/index/tximeta/linked_txome.json", **config),
            # symlink=expand(f"{{genome_dir}}/{{{{assembly}}}}/index/tximeta/{config['tximeta']['organism']}.{{{{assembly}}}}.{config['tximeta']['release']}.gtf", **config)
        params:
            source=config['tximeta']['source'],
            organism=config['tximeta']['organism'],
            release=config['tximeta']['release']
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-linked_txome.log", **config)
        conda:
            "../envs/gene_counts.yaml"
        resources:
            R_scripts=1 # conda's R can have issues when starting multiple times
        script:
            f"{config['rule_dir']}/../scripts/linked_txome.R"

    rule txi_count_matrix:
        """
        Convert estimated transcript abundances to gene count estimations and merge gene counts per assembly
        
        Also outputs a single cell experiment object similar to ARMOR (https://github.com/csoneson/ARMOR)
        
        Only works with Ensembl assemblies
        """
        input:
            linked_txome = expand("{genome_dir}/{{assembly}}/index/tximeta/linked_txome.json", **config),
            cts = get_counts
        output:
            counts = expand("{counts_dir}/{{assembly}}-counts.tsv", **config),
            SCE = expand("{counts_dir}/{{assembly}}-se.rds", **config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-txi_counts_matrix.log", **config)
        conda:
            "../envs/gene_counts.yaml"
        resources:
            R_scripts=1 # conda's R can have issues when starting multiple times
        script:
            f"{config['rule_dir']}/../scripts/txi.R"


elif config['quantifier'] == 'star':
    rule counts_matrix:
        """
        Merge gene counts per assembly
        
        From https://github.com/snakemake-workflows/rna-seq-star-deseq2
        """
        input:
            cts = get_counts
        output:
            expand("{counts_dir}/{{assembly}}-counts.tsv", **config)
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-counts_matrix.log", **config)
        run:
            import pandas as pd
            import os.path
            import sys

            with open(log[0], "w") as log_file:
                sys.stderr = sys.stdout = log_file

                # check if samples.tsv has at least one sample with strandedness
                if 'strandedness' in samples.columns:
                    def get_column(strandedness):
                        if pd.isnull(strandedness) or strandedness in ['no', 'NaN']:
                            return 1 #non stranded protocol
                        elif strandedness in ['forward', 'yes']:
                            return 2 #3rd column
                        elif strandedness == 'reverse':
                            return 3 #4th column, usually for Illumina truseq
                        else:
                            raise ValueError("'strandedness' column should be empty or have the " 
                                             "value 'none', 'yes', 'forward' (same as yes) or 'reverse', instead has the " 
                                             f"value {repr(strandedness)}")

                    # strandedness per sample/replicate
                    s2 = samples['strandedness']
                    if 'replicate' in samples:
                        s2 = samples.reset_index()[['replicate', 'strandedness']].drop_duplicates().set_index('replicate')

                counts = pd.DataFrame()
                for sample in input.cts:
                    sample_name = os.path.basename(sample).replace(wildcards.assembly + '-', '', 1)
                    strand_dependent_column = 1
                    if 'strandedness' in samples:
                        strandedness = s2["strandedness"].loc[sample_name]
                        strand_dependent_column = get_column(strandedness)

                    col = pd.read_csv(sample + '/ReadsPerGene.out.tab', sep='\t', index_col=0,
                                      usecols=[0, strand_dependent_column], header=None, skiprows=4)
                    col.columns = [sample_name]
                    counts = pd.concat([counts, col], axis = 1)

                counts.index.name = "gene"
                counts.to_csv(output[0], sep="\t")
