from io import StringIO
from contextlib import redirect_stdout, redirect_stderr

import pandas as pd

from seq2science.util import dense_samples


with open(str(snakemake.log), "w") as f:
    with redirect_stdout(f), redirect_stderr(f):
        counts = pd.read_table(snakemake.input[0], comment="#", index_col=0)

        samples = pd.read_table(StringIO(snakemake.params.samples), sep="\s+")
        samples.index.name = "sample"
        samples = samples.reset_index()
        samples = dense_samples(samples,
                                snakemake.params.technical_replicates,
                                snakemake.params.biological_replicates,
                                True,
                                False)
        breps = snakemake.params.breps

        groups = [breps.index(samples[samples["descriptive_name"] == trep]["biological_replicates"].values) for trep in counts.columns]

        newcounts = counts.groupby(groups, axis=1).mean()
        newcounts.columns = breps

        with open(snakemake.output[0], "w") as f:
            f.write(newcounts.to_csv(header=True, sep="\t"))
