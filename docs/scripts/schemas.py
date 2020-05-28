"""
This script auto-generates the
"""
import yaml
import os

final_md = (
"""
# All configurable options

TODO
This is an automatically generated ...
These are all configurable options, loosely spread around topics ....
"""
)

path = "schemas/config/"
order = {"general": "general",
         "alignment sth": "alignment_specific",
         "alignment": "alignment_general",
         "peak calling (ChIP & ATAC)": "peakcalling",
         "gene expression (RNA-seq)": "gene_expression",
         "trackhub": "trackhub"}

# we add 4 indentation as a start
indentation = -4

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
        final_md += ""
        final_md += f"{config}:\n"
        for key, val in settings.items():
            final_md = unpack_config(final_md, key, val, indentation)
        final_md += "```\n"
    final_md += "\n"

with open("docs/content/schemas.md", "w") as text_file:
    text_file.write(final_md)
