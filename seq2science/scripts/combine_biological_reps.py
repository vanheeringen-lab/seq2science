import pandas as pd

counts = pd.read_table(snakemake.input[0], comment="#", index_col=0)
print(counts)
local_samples = samples[samples["assembly"] == wildcards.assembly]
print(local_samples)

if descriptive_name in local_samples.columns:
   groups = [breps.index(local_samples[local_samples["descriptive_name"] == trep]["biological_replicates"].values[0]) for trep in counts.columns]
else:
   groups = [breps.index(local_samples[local_samples.index == trep]["biological_replicates"].values[0]) for trep in counts.columns]
print(groups)

newcounts = counts.groupby(groups, axis=1).mean()
print(newcounts)
newcounts.columns = breps
print(newcounts)

with open(snakemake.output[0]) as f:
   f.write(newcounts.to_csv(header=True, sep="\t"))
