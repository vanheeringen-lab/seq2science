"""
Make a pretty samples table and config yaml for in the MultiQC report
"""
from io import StringIO
import pandas as pd
from pretty_html_table import build_table


outstring = \
    "<!--\n" \
    "id: 'samplesconfig'\n" \
    "section_name: 'Samples & Config'\n" \
    "-->\n"

samples = pd.read_table(StringIO(snakemake.params.samples), sep="\s+")
outstring += "The samples file used for this run: <br>" \
             f"{build_table(samples, 'blue_dark')}"

if snakemake.params.config_used:
    outstring += "The config file used for this run: <br>"
    outstring += '<pre><code class="codeblock">'
    with open(snakemake.params.configfile, "r") as config_file:
        outstring += config_file.read()
    outstring += '</code></pre>'

with open(snakemake.output[0], "w") as out_file:
    out_file.write(outstring)
