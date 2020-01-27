def get_counts(wildcards):
    iterator = samples[samples['assembly'] == wildcards.assembly].index
    if config.get('combine_replicates', '') == 'merge' and 'condition' in samples:
        iterator = set(samples[samples['assembly'] == wildcards.assembly].condition)
    output = []
    for sample in iterator:
        output.append(f"{{result_dir}}/{{quantifier}}/{wildcards.assembly}-{sample}")
    return expand(output, **config)

if config['quantifier'] == 'salmon':
    def get_index(wildcards):
        index="{genome_dir}/" + wildcards.assembly + "/index/{quantifier}_decoy_aware" if config["decoy_aware_index"] else "{genome_dir}/" + wildcards.assembly + "/index/{quantifier}"
        return expand(index, **config)

    """
    currently does not work. Requires an updated version of tximeta on conda cloud ~1.4.4
    """
    rule linked_txome:
        """
        Generate a linked transcriptome for tximeta
        
        Also creates a symlink to the gtf in an Ensembl format (required by tximeta)
        
        Required to converting salmon output (estimated transcript abundances) to gene counts
        """
        input:
            fasta=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.gtf", **config),
            index_dir=get_index
        output:
            index=expand("{genome_dir}/{{assembly}}/index/tximeta/linked_txome.json", **config),
            symlink=expand(f"{{genome_dir}}/{{{{assembly}}}}/index/tximeta/{config['tximeta']['organism']}.{{{{assembly}}}}.{config['tximeta']['release']}.gtf", **config)
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
            "../scripts/linked_txome.R"

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
            "../scripts/txi.R"


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

            with open(log[0], "w") as log_file:
                def get_column(strandedness):
                    if pd.isnull(strandedness) or strandedness == 'no':
                        return 1 #non stranded protocol
                    elif strandedness in ['forward', 'yes']:
                        return 2 #3rd column
                    elif strandedness == 'reverse':
                        return 3 #4th column, usually for Illumina truseq
                    else:
                        raise ValueError(("'strandedness' column should be empty or have the " 
                                          "value 'none', 'yes', 'forward' (yes=forward) or 'reverse', instead has the " 
                                          "value {}").format(repr(strandedness)))

                counts = pd.DataFrame()
                for sample in input.cts:
                    sample_name = os.path.basename(sample).replace(wildcards.assembly + '-', '', 1)

                    # retrieve strandedness
                    if config.get('combine_replicates', '') == 'merge' and 'condition' in samples:
                        s2 = samples[['condition', 'strandedness']].drop_duplicates().set_index('condition')
                        strandedness = s2["strandedness"].loc[sample_name] if 'strandedness' in samples else 'no'
                    else:
                        strandedness = samples["strandedness"].loc[sample_name] if 'strandedness' in samples else 'no'

                    col = pd.read_csv(sample + '/ReadsPerGene.out.tab', sep='\t', index_col=0,
                                      usecols=[0, get_column(strandedness)], header=None, skiprows=4)
                    col.columns = [sample_name]
                    counts = pd.concat([counts, col], axis = 1)

                counts.index.name = "gene"
                counts.to_csv(output[0], sep="\t")
#
#
# else:
#     def get_counts(wildcards):
#         iterator = samples[samples['assembly'] == wildcards.assembly].index
#         if config.get('combine_replicates', '') == 'merge' and 'condition' in samples:
#             iterator = set(samples[samples['assembly'] == wildcards.assembly].condition)
#         output = []
#         for sample in iterator:
#             output.append(f"{{result_dir}}/{{quantifier}}/{wildcards.assembly}-{sample}.samtools-coordinate.bam")
#         return expand(output, **config)
#
#     rule counts_matrix:
#         input:
#             cts=get_counts,
#             gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.gtf", **config),
#         output:
#             expand("{result_dir}/gene_counts/{{assembly}}-counts.tsv", **config)
#         log:
#             expand("{log_dir}/counts_matrix/{{assembly}}-counts_matrix.log", **config)
#         conda:
#             "../envs/gene_counts.yaml"
#         shell:
#              """
#              htseq-count -f bam -r pos -m union {input.cts} {input.gtf} &> {log}
#              grep -v -- 'SAM\|GFF\|Warning' {log} > {output}
#              """
