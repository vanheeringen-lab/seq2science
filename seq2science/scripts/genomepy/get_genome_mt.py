import genomepy
import contextlib


annotation_path=snakemake.input[0]
logfile = snakemake.log[0]
output = snakemake.output[0]


def mito_genes(annotation_path):
    ann = genomepy.annotation.Annotation(annotation_path)
    gtf = ann.named_gtf
    
    # very quick, very dirty way to find the name for the mitochondrion.
    # could also try to read this from the assembly report, but that's also not perfect. 
    mt = gtf[gtf["seqname"].str.contains("chrM", case=False, regex=False)]["seqname"].unique()
    if len(mt) != 1:
        mt = gtf[gtf["seqname"].str.contains("MT", case=False, regex=False)]["seqname"].unique()
    if len(mt) != 1:
        mt = gtf[gtf["seqname"].str.contains("mito", case=False, regex=False)]["seqname"].unique()
    if len(mt) != 1:
        mt = gtf[gtf["seqname"].str.contains("m", case=False, regex=False)]["seqname"].unique()
    if len(mt) != 1:
        print("Could not extract mitochondrial gene list from annotation. We tried.....")
        return {'NONE'}
    
    genes = set(gtf[gtf["seqname"] == mt[0]].index)
    return genes


# redirect all messages to a logfile
with open(logfile, "w") as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        genomepy.logger.remove()
        genomepy.logger.add(
            logfile,
            format="<green>{time:HH:mm:ss}</green> <bold>|</bold> <blue>{level}</blue> <bold>|</bold> {message}",
            level="INFO",
        )
        try:
            genes = mito_genes(annotation_path)
            #Write genes to file
            with open(output, 'w') as f:
                for line in genes:
                    f.write(f"{line}\n")
        except Exception as e:
            print(e)
            print("\nSomething went wrong while extracting the MT genes from your annotation (see error message above). ")
