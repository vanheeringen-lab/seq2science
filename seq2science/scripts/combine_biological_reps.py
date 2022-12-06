from io import StringIO
from contextlib import redirect_stdout, redirect_stderr

import pandas as pd


with open(str(snakemake.log), "w") as f:
    with redirect_stdout(f), redirect_stderr(f):
        counts = pd.read_table(snakemake.input[0], comment="#", index_col=0)
        samples = pd.read_table(StringIO(snakemake.params.samples), sep="\s+")

        if "biological_replicates" in samples.columns:
            # set the index to the column with the names used in the counts table
            for col in ["descriptive_name", "technical_replicates"]:
                if col in samples.columns:
                    samples.set_index(col, inplace=True)
                    break

            # for each each column, determine which biological replicate it belongs to
            breps = snakemake.params.breps
            groups = []
            for column in counts.columns:
                idx = breps.index(samples.at[column, "biological_replicates"])
                groups.append(idx)

            # average all columns per biological replicate and assign the brep names to the results
            newcounts = counts.groupby(groups, axis=1).mean()
            newcounts.columns = breps

        else:
            newcounts = counts.copy()

        with open(snakemake.output[0], "w") as f:
            f.write(newcounts.to_csv(header=True, sep="\t"))
