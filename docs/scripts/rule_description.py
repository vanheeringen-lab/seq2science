"""
This script auto-generates the
"""
import yaml
import os
import re

final_md = (
"""
# Per rule explanation

This is an automatically generated list of all supported rules, and their docstrings. At the start of each workflow 
run a list is printed of which rules will be run. And while the workflow is running it prints which rules are being \
started and finished. This page is here to give an explanation to the user about what each rule does, and for \
developers to find what is, and isn't yet supported.
"""
)

path = "rules/"


def get_dirty_docstrings(string):
    splitter = re.compile("rule (.*):[\s\S]*?\"\"\"([\s\S]*?)\"\"\"", re.MULTILINE)
    docstrings = {}
    for match in splitter.finditer(string):
        docstrings[match.group(1)] = match.group(2)
    return docstrings


def get_clean_docstrings(dirty):
    clean = {}
    for rule, docstring in dirty.items():
        firstline = docstring.split("\n")[1]

        indentation = len(firstline) - len(firstline.lstrip())
        docstring = docstring.replace(" " * indentation, "")
        docstring = docstring.replace(" " * (indentation - 4), "")
        docstring = docstring.strip("\n")
        clean[rule] = docstring

    return clean

all_rules = {}
for rules_file in os.listdir(path):
    with open(path + rules_file, 'r') as file:
        dirty = get_dirty_docstrings(file.read())
    clean = get_clean_docstrings(dirty)
    all_rules.update(clean)

for rule in sorted(all_rules.keys()):
    docstring = all_rules[rule]

    final_md += f"####{rule}\n"
    final_md += "```\n"
    final_md += f"{docstring}\n"
    final_md += "```\n"
    final_md += f"\n"

with open("docs/content/all_rules.md", "w") as text_file:
    text_file.write(final_md)
