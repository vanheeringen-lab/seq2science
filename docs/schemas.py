import yaml
import os

final_md = (
"""
# All options

This is an automatically generated ...
"""
)

path = "schemas/config/"
order = {"general": "general",
         "alignment sth": "alignment_specific",
         "alignment": "alignment_general",
         "peak calling (ChIP & ATAC)": "peakcalling",
         "gene expression (RNA-seq)": "gene_expression",
         "trackhub": "trackhub"}

for name, file in order.items():
    with open(path + file + ".schema.yaml", 'r') as stream:
        schema = yaml.safe_load(stream)
    final_md += f"## {name}\n"
    for config, settings in schema["properties"].items():
        final_md += "```\n"
        final_md += f"{config}:\n"
        for key, val in settings.items():
            final_md += f"{key:>15}: {val:}\n"
        final_md += "```\n"
    final_md += "\n"

with open("docs/schemas.md", "w") as text_file:
    text_file.write(final_md)
