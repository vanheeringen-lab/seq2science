import sys

blacklisted_rule_names = [
    "multiqc_header_info",
    "multiqc_rename_buttons",
    "multiqc_schema",
    "multiqc_samplesconfig",
    "multiqc_explain",
    "multiqc_filter_buttons",
]

graph_file = sys.argv[1]
with open(graph_file) as f:
    lines = f.readlines()

blacklisted_node_ids = []
for line in lines:
    for rule in blacklisted_rule_names:
        if rule in line:
            node_id = line.strip().split("[")[0]
            blacklisted_node_ids.append(node_id)

with open(graph_file, "w") as f:
    for line in lines:
        line_start = line.strip().split()[0]
        for node_id in blacklisted_node_ids:
            if line_start.startswith(node_id):
                # line contains a blacklisted node ID.
                break
        else:
            f.write(line)
