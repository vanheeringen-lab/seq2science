def get_counts(wildcards):
    output = []
    for sample in samples[samples['assembly'] == wildcards.assembly].index:
        output.append(f"{{result_dir}}/{{aligner}}/{wildcards.assembly}-{sample}")
    return expand(output, **config)

if config['aligner'] == 'salmon':
    rule linked_txome:
        """
        Generate a linked transcriptome for tximeta
        
        Also creates a symlink to the gtf in an Ensembl format (required by tximeta)
          
        Required to converting salmon output (estimated transcript abundances) to gene counts
        """
        input:
            fasta=expand("{genome_dir}/{{assembly}}/{{assembly}}.transcripts.fa", **config),
            gtf=expand("{genome_dir}/{{assembly}}/{{assembly}}.gtf", **config),
            index_dir=expand("{genome_dir}/{{assembly}}/index/{aligner}", **config)
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
            "../envs/txi.yaml"
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
            counts = expand("{result_dir}/gene_counts/{{assembly}}-counts.tsv", **config),
            SCE = expand("{result_dir}/gene_counts/{{assembly}}-se.rds", **config),
        log:
            expand("{log_dir}/counts_matrix/{{assembly}}-txi_counts_matrix.log", **config)
        conda:
            "../envs/txi.yaml"
        script:
            "../scripts/txi.R"


elif config['aligner'] == 'star':
    rule counts_matrix:
        """
        Merge gene counts per assembly
        
        From https://github.com/snakemake-workflows/rna-seq-star-deseq2
        """
        input:
            cts = get_counts
        output:
            expand("{result_dir}/gene_counts/{{assembly}}-counts.tsv", **config)
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
                    strandedness = samples["strandedness"].loc[sample_name] if 'strandedness' in samples else 'no'
                    col = pd.read_csv(sample + '/ReadsPerGene.out.tab', sep='\t', index_col=0,
                                      usecols=[0, get_column(strandedness)], header=None, skiprows=4)
                    col.columns = [sample_name]
                    counts = pd.concat([counts, col], axis = 1)

                counts.index.name = "gene"
                counts.to_csv(output[0], sep="\t")
