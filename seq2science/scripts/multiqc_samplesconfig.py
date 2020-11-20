"""
Make a pretty samples table and config yaml for in the MultiQC report
"""
from pretty_html_table import build_table


outstring = \
    "<!--\n" \
    "id: 'samplesconfig'\n" \
    "section_name: 'Samples & Config'\n" \
    "-->\n"


outstring += "The samples file used for this run: <br>" \
             f"{build_table(snakemake.params.sanitized_samples, 'blue_dark')}"

if snakemake.params.config_used:
    outstring += "The config file used for this run: <br>"
    outstring += '<pre><code class="codeblock">'
    with open(snakemake.params.configfile, "r") as config_file:
        outstring += config_file.read()
    outstring += '</code></pre>'

with open(snakemake.output[0], "w") as out_file:
    out_file.write(outstring)
