import re
import networkx as nx

rules_to_keep = {
    # "0.13 0.6 0.85",  # yellow
    # "0.09 0.6 0.85",  # brown
    # "0.28 0.6 0.85",  # green
    # "0.04 0.6 0.85",  # red
    # "0.00 0.6 0.85",  # cherry
    # "0.63 0.6 0.85",  # purple
    # "0.59 0.6 0.85",  # dark blue
    # "0.58 0.6 0.85",  # blue
    # "0.49 0.6 0.85",  # teal
    # input
    "get_genome": "0.49 0.6 0.85",  # teal
    "ena2fastq_SE": "0.49 0.6 0.85",  # teal
    "ena2fastq_PE": "0.49 0.6 0.85",  # teal
    "sra2fastq_SE": "0.49 0.6 0.85",  # teal
    "sra2fastq_PE": "0.49 0.6 0.85",  # teal
    # fastq
    "fastp_SE": "0.13 0.6 0.85",  # yellow
    "fastp_PE": "0.13 0.6 0.85",  # yellow
    "trimgalore_SE": "0.13 0.6 0.85",  # yellow
    "trimgalore_PE": "0.13 0.6 0.85",  # yellow
    "merge_replicates": "0.13 0.6 0.85",  # yellow
    # align
    "bowtie2_align": "0.13 0.6 0.85",  # yellow
    "bwa_mem": "0.13 0.6 0.85",  # yellow
    "bwa_mem2": "0.13 0.6 0.85",  # yellow
    "hisat2_align": "0.13 0.6 0.85",  # yellow
    "minimap2_align": "0.13 0.6 0.85",  # yellow
    "star_align": "0.13 0.6 0.85",  # yellow
    "mark_duplicates": "0.13 0.6 0.85",  # yellow
    "sieve_bam": "0.13 0.6 0.85",  # yellow
    # peak counting
    "macs2_callpeak": "0.13 0.6 0.85",  # yellow
    "call_peak_genrich": "0.13 0.6 0.85",  # yellow
    "hmmratac": "0.13 0.6 0.85",  # yellow
    "create_SNAP_object": "0.13 0.6 0.85",  # yellow
    # gene counting/quantification
    "htseq_count": "0.13 0.6 0.85",  # yellow
    "featurecounts": "0.13 0.6 0.85",  # yellow
    "salmon_quant": "0.13 0.6 0.85",  # yellow
    # trackhub
    "bam_bigwig": "0.00 0.6 0.85",  # cherry
    "peak_bigpeak": "0.00 0.6 0.85",  # cherry
    "bedgraph_bigwig": "0.00 0.6 0.85",  # cherry
    "trackhub": "0.00 0.6 0.85",  # cherry
    # multiqc
    "multiqc": "0.63 0.6 0.85",  # purple
    # peak files
    "coverage_table": "0.28 0.6 0.85",  # green
    "onehot_peaks": "0.28 0.6 0.85",  # green
    "create_bins_SNAP_object": "0.28 0.6 0.85",  # green
    # gene files
    "gene_id2name": "0.28 0.6 0.85",  # green
    "tpm_matrix": "0.28 0.6 0.85",  # green
    "count_matrix": "0.28 0.6 0.85",  # green
    "txi_count_matrix": "0.28 0.6 0.85",  # green
    "pytxi_count_matrix": "0.28 0.6 0.85",  # green
    "citeseqcount": "0.28 0.6 0.85",  # green
    "kallistobus_count": "0.28 0.6 0.85",  # green
    # other
    "gimme_maelstrom": "0.28 0.6 0.85",  # green
    "deseq2": "0.28 0.6 0.85",  # green
    "dexseq_count_matrix": "0.28 0.6 0.85",  # green
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
        for parent in parents:
            for daughter in daughters:
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
