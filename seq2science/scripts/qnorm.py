import pandas as pd
import qnorm


df = pd.read_csv(snakemake.input[0], comment="#", index_col=0, sep="\t")

# cpm normalization
df = df * 1_000_000 / df.sum(axis=0)

# quantile normalize
df_qn = qnorm.quantile_normalize(df)
open(str(snakemake.output[0]), "w").write(
    "# The number of reads under each peak, cpm quantile normalized\n" +
    df_qn.to_csv(index_label="loc", index=True, header=True, sep="\t")
)
