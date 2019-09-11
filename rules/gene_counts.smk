def get_counts(wildcards):
    output = []
    for fname in samples[samples['assembly'] == wildcards.assembly].index:
        output.append(f"{{result_dir}}/{{aligner}}/{wildcards.assembly}-{fname}")
    return expand(output, **config)

rule count_matrix:
    """
    Merge gene counts per assembly
    
    From https://github.com/snakemake-workflows/rna-seq-star-deseq2
    """
    input:
        cts = get_counts
    output:
        expand("{result_dir}/gene_counts/{{assembly}}-counts.tsv", **config)
    run:
        import pandas as pd
        import os.path

        def get_column(strandedness):
            if pd.isnull(strandedness) or strandedness in ["none", "no", "unknown"]:
                return 1 #non stranded protocol
            elif strandedness in ["forward", 'yes']:
                return 2 #3rd column
            elif strandedness == "reverse":
                return 3 #4th column, usually for Illumina truseq
            else:
                raise ValueError(("'strandedness' column should be empty or have the " 
                                  "value 'none', 'yes', 'forward' (yes=forward) or 'reverse', instead has the " 
                                  "value {}").format(repr(strandedness)))

        counts = pd.DataFrame()
        for sample in input.cts:
            sample_name = os.path.basename(sample).replace(wildcards.assembly + '-', '', 1)
            strandedness = samples["strandedness"].loc[sample_name]
            col = pd.read_csv(sample + '/ReadsPerGene.out.tab', sep='\t', index_col=0,
                              usecols=[0, get_column(strandedness)], header=None, skiprows=4)
            col.columns = [sample_name]
            counts = pd.concat([counts, col], axis = 1)

        counts.index.name = "gene"
        counts.to_csv(output[0], sep="\t")
