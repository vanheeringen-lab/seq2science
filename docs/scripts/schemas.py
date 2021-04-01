"""
This script auto-generates the
"""
import yaml
import os

final_md = (
"""
# All configurable options

This is an automatically generated summary of all configurable options for seq2science. These options are loosely \
grouped around workflows/topics, however they are generally also shared across workflows. So it is possible that  \
tunable configuration settings are not mentioned in their topic. At the start of each seq2science run the complete \
configuration is printed to stdout. You can use that printed configuration as the complete list of tunable \
configuration settings. 

We believe that all our default settings are reasonable, and manual finetuning is generally not required.
"""
)

path = "seq2science/schemas/config/"
order = {
    "General": "general",
    "Download": "download",
    "Alignment general": "alignment_general",
    "Workflow: Alignment": "alignment_specific",
    "Workflow: ChIP & ATAC-seq": "peakcalling",
    "Workflow: RNA-seq": "gene_expression",
    "Workflow: Single-cell ATAC-seq": "scatac",
    "Workflow: Single-cell RNA-seq": "scrna",
    "Differential gene/peak analysis": "deseq2",
    "Trackhub": "trackhub",
}

# we add 4 indentation as a start
indentation = 0

def unpack_config(markdown, key, val, indentation):
    indentation += 4
    spaces = ' ' * indentation
    if isinstance(val, dict):
        markdown += f"{spaces}{key}:\n"
        for newkey, newval in val.items():
            markdown = unpack_config(markdown, newkey, newval, indentation)
    else:
        markdown += f"{spaces}{key}: {val}\n"
    return markdown

for name, file in order.items():
    with open(path + file + ".schema.yaml", 'r') as stream:
        schema = yaml.safe_load(stream)
    final_md += f"## {name}\n"
    for config, settings in schema["properties"].items():
        final_md += f"#####{config}\n```\n"
        final_md += f"{config}:\n"
        for key, val in settings.items():
            final_md = unpack_config(final_md, key, val, indentation)
        final_md += "```\n"
    final_md += "\n"

with open("docs/content/schemas.md", "w") as text_file:
    text_file.write(final_md)
