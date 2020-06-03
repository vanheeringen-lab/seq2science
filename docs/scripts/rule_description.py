"""
This script auto-generates the
"""
import yaml
import os
import re

final_md = (
"""\
# Per rule explanation

This is an automatically generated list of all supported rules, their docstrings, and command. At the start of each \
workflow run a list is printed of which rules will be run. And while the workflow is running it prints which rules are \
being started and finished. This page is here to give an explanation to the user about what each rule does, and for \
developers to find what is, and isn't yet supported.

"""
)

path = "seq2science/rules/"


def get_dirty_docstrings(string):
    splitter = re.compile("rule (.*):[\s\S]*?\"\"\"([\s\S]*?)\"\"\"", re.MULTILINE)
    docstrings = {}
    for match in splitter.finditer(string):
        docstrings[match.group(1)] = match.group(2)
    return docstrings


def cleanup(dirty):
    clean = {}
    for rule, docstring in dirty.items():
        firstline = docstring.split("\n")[1]

        indentation = len(firstline) - len(firstline.lstrip())
        docstring = docstring.replace(" " * indentation, "")
        docstring = docstring.replace(" " * (indentation - 4), "")
        docstring = docstring.strip("\n")
        clean[rule] = docstring

    return clean


def get_dirty_shell(string):
    splitter = re.compile("rule (.*):[\s\S]*?shell:[\s\S]*?\"\"\"[\s\S]([\s\S]*?)\"\"\"", re.MULTILINE)
    shell_cmds = {}
    for match in splitter.finditer(string):
        shell_cmds[match.group(1)] = match.group(2)
    return shell_cmds


all_rules_doc = {}
all_rules_shell = {}
for rules_file in os.listdir(path):
    with open(path + rules_file, 'r') as file:
        text = file.read()
    shell_cmd = cleanup(get_dirty_shell(text))
    all_rules_shell.update(shell_cmd)

    docstrings = cleanup(get_dirty_docstrings(text))
    all_rules_doc.update(docstrings)

for rule in sorted(all_rules_doc.keys()):
    docstring = all_rules_doc[rule]

    final_md += f"#### {rule}\n"
    final_md += f"{docstring}\n"
    if rule in all_rules_shell:
        final_md += "```\n"
        final_md += f"{all_rules_shell[rule]}\n"
        final_md += "```\n"
    final_md += f"\n"

with open("docs/content/all_rules.md", "w") as text_file:
    text_file.write(final_md)
