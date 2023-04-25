import re
import networkx as nx


colors = {
    "yellow": "0.13 0.6 0.85",
    "brown": "0.09 0.6 0.85",
    "green": "0.28 0.6 0.85",
    "red": "0.04 0.6 0.85",
    "cherry": "0.00 0.6 0.85",
    "purple": "0.63 0.6 0.85",
    "dark blue": "0.59 0.6 0.85",
    "blue": "0.58 0.6 0.85",
    "teal": "0.49 0.6 0.85",
}

rules_to_keep = {
    # input
    "get_genome": colors["teal"],
    "ena2fastq_SE": colors["teal"],
    "ena2fastq_PE": colors["teal"],
    "sra2fastq_SE": colors["teal"],
    "sra2fastq_PE": colors["teal"],

    # fastq
    "fastp_SE": colors["yellow"],
    "fastp_PE": colors["yellow"],
    "trimgalore_SE": colors["yellow"],
    "trimgalore_PE": colors["yellow"],
    # "merge_replicates": colors["yellow"],

    # alignment
    "bowtie2_align": colors["yellow"],
    "bwa_mem": colors["yellow"],
    "bwa_mem2": colors["yellow"],
    "hisat2_align": colors["yellow"],
    "minimap2_align": colors["yellow"],
    "star_align": colors["yellow"],
    "mark_duplicates": colors["yellow"],
    # "sieve_bam": colors["yellow"],

    # peak counting
    "macs2_callpeak": colors["yellow"],
    "call_peak_genrich": colors["yellow"],
    "hmmratac": colors["yellow"],
    "create_SNAP_object": colors["yellow"],

    # gene counting/quantification
    "htseq_count": colors["yellow"],
    "featurecounts": colors["yellow"],
    "salmon_quant": colors["yellow"],

    # trackhub
    # "bam_bigwig": colors["cherry"],
    # "peak_bigpeak": colors["cherry"],
    # "bedgraph_bigwig": colors["cherry"],
    "trackhub": colors["cherry"],

    # multiqc
    "multiqc": colors["purple"],

    # peak files
    "coverage_table": colors["green"],
    # "onehot_peaks": colors["green"],
    "create_bins_SNAP_object": colors["green"],

    # gene files
    # "gene_id2name": colors["green"],
    "tpm_matrix": colors["green"],
    "count_matrix": colors["green"],
    "txi_count_matrix": colors["green"],
    "pytxi_count_matrix": colors["green"],
    "citeseqcount": colors["green"],
    "kallistobus_count": colors["green"],

    # analysis
    "gimme_maelstrom": colors["green"],
    "deseq2": colors["green"],
    "dexseq_count_matrix": colors["green"],
}


class Digraph:
    def __init__(self, infile):
        with open(infile) as f:
            lines = f.readlines()
        self.type, self.name = lines[0].split()[0:2]
        self.graph_style = lines[1]
        self.node_style = lines[2]
        self.edge_style = lines[3]

        self.nodes = dict()
        self.edges = set()
        self.label2id = dict()
        l = re.compile(r'label = "(.*?)"')
        c = re.compile(r'color = "(.*?)"')
        s = re.compile(r'style="(.*?)"')
        for line in lines[4:]:
            line = line.strip()

            # read edges
            edge = tuple(line.split(" -> "))
            if len(edge) == 2:
                self.edges.add(edge)

            # read nodes
            elif "[" in line[:5]:
                node_id = line.split("[")[0]
                label = l.search(line).groups()[0]
                color = c.search(line).groups()[0]
                style = s.search(line).groups()[0]
                self.nodes[node_id] = {
                    "label": label,
                    "color": color,
                    "style": style,
                }
                self.label2id[label] = node_id

    def _order_edges(self):
        """
        edges are sorted by
        1) ascending target nodes and
        2) descending source nodes
        """
        ordered = []
        sources = sorted(set(int(e[0]) for e in list(self.edges)), reverse=True)
        targets = sorted(set(int(e[1]) for e in list(self.edges)))
        for target in targets:
            for source in sources:
                edge = (str(source), str(target))
                if edge in self.edges:
                    ordered.append(edge)
        return ordered

    def write(self, fname):
        with open(fname, "w") as f:
            f.write(" ".join([self.type, self.name, "{\n"]))
            f.write(self.graph_style)
            f.write(self.node_style)
            f.write(self.edge_style)
            for k, v in self.nodes.items():
                l = v["label"]
                c = v["color"]
                s = v["style"]
                line = f'    {k}[label = "{l}", color = "{c}", style="{s}"];\n'
                f.write(line)
            for a, b in self._order_edges():
                line = f"    {a} -> {b}\n"
                f.write(line)
            f.write("}\n")

    def _get_node_id(self, node):
        node = str(node)
        if node.isdigit() and node in self.nodes:
            return node
        return self.label2id.get(node)

    def color_node(self, node, color):
        node_id = self._get_node_id(node)
        if node_id is None:
            return

        self.nodes[node_id]["color"] = color

    def remove_node(self, node):
        node_id = self._get_node_id(node)
        if node_id is None:
            return

        # remove the node
        label = self.nodes[node_id]["label"]
        del self.nodes[node_id]
        del self.label2id[label]
        # remove all edges with this node
        for edge in self.edges.copy():
            if node_id in edge:
                self.edges.remove(edge)

    def remove_edge(self, source, target):
        source_id = self._get_node_id(source)
        target_id = self._get_node_id(target)
        edge = (source_id, target_id)
        if edge in self.edges:
            self.edges.remove(edge)

    def hide_node(self, node):
        """remove a node, and connect incoming edges and outgoing edges"""
        node_id = self._get_node_id(node)
        if node_id is None:
            return

        # identify parent and daughter nodes
        parents = []
        daughters = []
        for edge in self.edges:
            if node_id == edge[0]:
                daughters.append(edge[1])
            elif node_id == edge[1]:
                parents.append(edge[0])

        # remove the node
        self.remove_node(node)
        # connect the neighboring nodes
        scrna = len({"citeseqcount", "kallistobus_count"} & {self.nodes[n]["label"] for n in self.nodes}) > 0
        for parent in parents:
            for daughter in daughters:
                # don't mix paired end and single end nodes (unless its scRNA-seq)
                p = self.nodes[parent]["label"][-3:]
                d = self.nodes[daughter]["label"][-3:]
                if (p == "_SE" and d == "_PE") or (p == "_PE" and d == "_SE") and not scrna:
                    continue
                edge = (parent, daughter)
                self.edges.add(edge)

    def transitive_reduction(self):
        g = nx.DiGraph(self.edges)
        g = nx.algorithms.transitive_reduction(g)
        self.edges = set(g.edges())

    # def _get_edges(self, node_id, kind="any"):
    #     if kind == "any":
    #         return [e for e in self.edges if node_id in e]
    #     if kind == "parents":
    #         return [e for e in self.edges if e[1] == node_id]
    #     if kind == "daughters":
    #         return [e for e in self.edges if e[0] == node_id]
    #     raise ValueError
