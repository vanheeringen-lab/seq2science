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
upsetplot = plot(from_contents(test_data), fig=f)

plt.savefig(snakemake.output[0], dpi=450)
