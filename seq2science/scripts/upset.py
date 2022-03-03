import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import from_contents, from_memberships, from_indicators, plot


# read the table
df = pd.read_table(snakemake.input[0], comment="#", index_col=0)
df = df.drop(":-")
df = df > 0

# put it in a dictionary
data = dict()
for col in df.columns:
    row = list()
    for index, val in df[col].iteritems():
        if val:
            row.append(index)
    data[col] = row


f, ax = plt.subplots()
ax.axis("off")

# generate the upset plot
# but make sure to filter to a max of 31 combinations
# 31 is the maximum of combinations of 5 different items
content = from_contents(data)
uniques, counts = np.unique(content.index, return_counts=True)

sorted_uniques = [x for _, x in sorted(zip(counts, uniques), reverse=True)]

plot(content.loc[sorted_uniques[:31]], sort_by=None)

upsetplot = plot(from_contents(data), fig=f)

plt.savefig(snakemake.output[0], dpi=250)
