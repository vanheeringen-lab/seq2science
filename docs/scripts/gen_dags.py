"""
This script auto-generates workflow DAGs in docs/resources/*.png
"""
import os
from os.path import abspath, join
from pathlib import Path
import subprocess as sp

from clean_dags import rules_to_hide, rules_to_color, Digraph


rule_dir = abspath("seq2science/rules")
out_dir = abspath("docs/resources")
in_dir = abspath("seq2science/workflows")
for workflow in os.listdir(in_dir):
    workflow_dir = join(in_dir, workflow)

    snakefile = join(workflow_dir, "Snakefile")
    samples = join(workflow_dir, "samples.tsv")
    config = join(workflow_dir, "config.yaml")
    tmp = join(workflow_dir, ".tmp_graph.txt")
    graph = join(out_dir, f"{workflow}.png")
    if workflow == "scrna_seq":
        barcodes = join(workflow_dir, 'barcodes.txt')
        Path(barcodes).touch()

    # create a DAG rulegraph
    sp.check_call(
        f'snakemake -s {snakefile} --configfile {config} --config samples={samples} rule_dir={rule_dir} '
        f'--quiet --rulegraph > {tmp}', shell=True
    )

    # clean the DAG rulegraph
    g = Digraph(tmp)
    for rule in rules_to_hide:
        g.hide_node(rule)
    for rule in g.nodes:
        g.nodes[rule]["color"] = "0.13 0.6 0.85"  # yellow
    for rule, color in rules_to_color.items():
        g.color_node(rule, color)
    g.transitive_reduction()
    g.write(tmp)

    # convert to image
    sp.check_call(f'dot -Tpng -Gbgcolor=transparent -Gdpi=450 {tmp} > {graph}', shell=True)

    os.remove(tmp)
