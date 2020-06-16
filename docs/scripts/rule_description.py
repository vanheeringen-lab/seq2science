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


def cleanup_docstring(dirty):
    clean = {}
    for rule, docstring in dirty.items():
        firstline = docstring.split("\n")[1]

        indentation = len(firstline) - len(firstline.lstrip())
        docstring = docstring.replace(" " * indentation, "")
        docstring = docstring.replace(" " * (indentation - 4), "")
        docstring = docstring.strip("\n")
        clean[rule] = docstring

    return clean


def cleanup_shell(dirty):
    clean = {}
    for rule, shell in dirty.items():
        firstline = shell.split("\n")[1]

        indentation = len(firstline) - len(firstline.lstrip())
        docstring = "\n".join([shell_line.replace(" " * indentation, "", 1) for shell_line in shell.split("\n")])
        docstring = docstring.strip("\n")
        clean[rule] = docstring

    return clean


def get_dirty_shell(string):
    splitter = re.compile("rule (.*):[\s\S]*?shell:[\s\S]*?\"\"\"[\s\S]([\s\S]*?)\"\"\"", re.MULTILINE)
    shell_cmds = {}
    for substring in string.split("\n\n\n"):
        for match in splitter.finditer(substring):
            shell_cmds[match.group(1)] = match.group(2)
    return shell_cmds


all_rules_doc = {}
all_rules_shell = {}
for rules_file in os.listdir(path):
    with open(path + rules_file, 'r') as file:
        text = file.read()
    shell_cmd = cleanup_shell(get_dirty_shell(text))
    all_rules_shell.update(shell_cmd)

    docstrings = cleanup_docstring(get_dirty_docstrings(text))
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
