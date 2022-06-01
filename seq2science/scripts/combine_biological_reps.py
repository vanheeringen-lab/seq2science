from io import StringIO
from contextlib import redirect_stdout, redirect_stderr

import pandas as pd

with open(str(snakemake.log), "w") as f:
   with redirect_stdout(f), redirect_stderr(f):
      counts = pd.read_table(snakemake.input[0], comment="#", index_col=0)
      samples = pd.read_table(StringIO(snakemake.params.samples), sep="\s+")
      breps = snakemake.params.breps

      print(breps)
      print(samples)
      print(counts.head())
      print(samples.columns)
      for trep in counts.columns:
          print(trep, samples[samples["descriptive_name"] == trep]["biological_replicates"].values)
      if "descriptive_name" in samples.columns:
         groups = [breps.index(samples[samples["descriptive_name"] == trep]["biological_replicates"].values) for trep in counts.columns]
      else:
         groups = [breps.index(samples[samples.index == trep]["biological_replicates"].values[0]) for trep in counts.columns]

      newcounts = counts.groupby(groups, axis=1).mean()
      newcounts.columns = breps

      with open(snakemake.output[0], "w") as f:
         f.write(newcounts.to_csv(header=True, sep="\t"))
