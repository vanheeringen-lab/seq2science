import io
import base64
from collections import Counter
from itertools import accumulate

import pyfaidx
import matplotlib.pyplot as plt


outfile = snakemake.output[0]
assembly = snakemake.wildcards.assembly

# empty old file (if exists)
open(outfile, "w").close()

# read the assembly
fa = pyfaidx.Fasta(snakemake.input.genome)

# variables to keep track of
gc = at = 0
sizes = list()

# count gc and at occurances and check how long each contig is
for key in fa.keys():
    sizes += [len(fa[key])]
    c = Counter(str(fa[key]))
    gc += c.get("G", 0) + c.get("C", 0) + c.get("g", 0) + c.get("c", 0)
    at += c.get("A", 0) + c.get("T", 0) + c.get("a", 0) + c.get("t", 0)

# now sort the sizes from large to small
sizes = list(sorted(sizes, reverse=True))

# get the cumulative and total sizes
cum_sizes = list(accumulate(sizes))
total_size = sum(sizes)

# get the n50 and l50 stats
n50_done = n75_done = False
for i, size in enumerate(cum_sizes):
    if not n50_done and size > total_size / 100 * 50:
        n50_done = True
        n50 = (i + 1, sizes[i])
    if not n75_done and size > total_size / 100 * 75:
        n75_done = True
        n75 = (i + 1, sizes[i])

# make an image off the size distribution of contigs
# make the actual plot
fig, ax = plt.subplots()
ax.plot(sizes)
ax.set_title("Contig size distribution")
ax.set_xlabel("contig number (ordered by size)")
ax.set_ylabel("contig size")
plt.xscale("log")

# now save the image as html text
img = io.BytesIO()
fig.savefig(img, format='png')
img.seek(0)
html = '<img src="data:image/png;base64, {}">'.format(base64.b64encode(img.getvalue()).decode('utf-8'))

# if we have an annotation check for the number of genes present
if hasattr(snakemake.input, "annotation"):
    with open(snakemake.input.annotation) as f:
        annotation = f.read()
    nr_genes = annotation.count("gene\t")
    annotation_text = f"""The genome annotation contains {nr_genes} genes."""
else:
    annotation_text = ""

# save it all
with open(outfile, "a") as f:
    f.write(
        f"""
<!--
id: 'assembly_stats'
section_name: 'Assembly stats'
-->
Genome assembly {assembly} contains of {len(sizes)} contigs, with a GC-to-AT ratio \
of {gc / (gc + at) * 100:.2f}%, and {(total_size - gc - at) / total_size * 100:.2f}%\
 consists of the letter N. The <a href="https://en.wikipedia.org/wiki/N50,_L50,_and_r\
elated_statistics">N50-L50</a> stats are {n50[1]}-{n50[0]} and the N75-L75 stats are \
{n75[1]}-{n75[0]}. {annotation_text}

<br>

<center>
{html}
</center>
"""
    )
