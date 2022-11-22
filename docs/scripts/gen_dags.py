"""
This script auto-generates workflow DAGs in docs/resources/*.png
"""
import os
from os.path import abspath, join
from pathlib import Path
import sys
import subprocess as sp

from clean_dags import rules_to_keep, Digraph


conda_dir = sys.base_exec_prefix
rule_dir = abspath("seq2science/rules")
out_dir = abspath("docs/resources")
in_dir = abspath("seq2science/workflows")
for workflow in os.listdir(in_dir):
    workflow_dir = join(in_dir, workflow)
    if workflow == "scrna_seq":
        barcodes = join(workflow_dir, "barcodes.txt")
        Path(barcodes).touch()

    # create a DAG rulegraph
    cmd = join(conda_dir, "bin", "snakemake")
    snakefile = join(workflow_dir, "Snakefile")
    samples = join(workflow_dir, "samples.tsv")
    config = join(workflow_dir, "config.yaml")
    tmp = join(workflow_dir, ".tmp_graph.txt")
    sp.check_call(
        f"{cmd} -s {snakefile} --configfile {config} "
        f"--config samples={samples} rule_dir={rule_dir} "
        f"--quiet --rulegraph > {tmp}",
        shell=True,
    )

    # clean the DAG rulegraph
    g = Digraph(tmp)
    for node in g.nodes.copy():
        rule = g.nodes[node]["label"]
        if rule in rules_to_keep:
            g.nodes[node]["color"] = rules_to_keep[rule]
        else:
            g.hide_node(rule)
    g.transitive_reduction()
    g.write(tmp)

    # convert to image
    cmd = join(conda_dir, "bin", "dot")
    graph = join(out_dir, f"{workflow}.png")
    sp.check_call(
        f"{cmd} -Tpng -Gbgcolor=transparent -Gdpi=450 {tmp} > {graph}", shell=True,
    )

    os.remove(tmp)
