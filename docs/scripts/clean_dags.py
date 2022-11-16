import re
import networkx as nx

rules_to_hide = [
    # preprocessing
    "runs2sample",
    "run2sra",
    "fastq_pair",
    "get_genome_support_files",
    "get_genome_annotation",
    "get_genome_blacklist",
    "setup_blacklist",
    "complement_blacklist",
    "extend_genome_annotation",
    "extend_genome_blacklist",
    "get_effective_genome_size",

    # aligner indexes
    "bowtie2_index",
    "bwa_index",
    "bwa_mem2_index",
    "hisat2_splice_aware_index",
    "hisat2_index",
    "minimap2_index",
    "star_index",

    # quantifier indexes
    "get_transcripts",
    "partial_decoy_transcripts",
    "full_decoy_transcripts",
    "salmon_index",
    "linked_txome",
    "kallistobus_ref",
    "kallistobus_ref_kite",

    # sam/bam intermediate rules
    "samtools_presort",
    "samtools_index",
    "samtools_sort",
    "sambamba_sort",
    "samtools_sort_allsizes",
    "sam2bam",

    # genes/peaks intermediate rules
    "macs_bdgcmp",
    "macs_cmbreps",
    "idr",
    "genrich_pileup",
    "keep_mates",
    "hmmratac_genome_info",
    "narrowpeak_summit",
    "combine_peaks",
    "random_subset_peaks",
    "bedtools_slop",
    "motif2factors",
    "prepare_DEXseq_annotation",
    "log_normalization",
    "quantile_normalization",
    "edgeR_normalization",
    "mean_center",
    "combine_biological_reps",
    # "onehot_peaks",

    # trackhub input
    "softmask_track_1",
    "softmask_track_2",
    "trackhub_index",
    "twobit",
    "gcPercent",
    "cytoband",

    # multiqc input
    "multiqc_header_info",
    "multiqc_rename_buttons",
    "multiqc_schema",
    "multiqc_samplesconfig",
    "multiqc_explain",
    "multiqc_filter_buttons",
    "multiqc_assembly_stats",
    "fastp_qc_SE",
    "fastp_qc_PE",
    "fastqc",
    "mt_nuc_ratio_calculator",
    "samtools_stats",
    "insert_size_metrics",
    "multiBamSummary",
    "plotPCA",
    "plotCorrelation",
    "plotFingerprint",
    "computeMatrix_peak",
    "computeMatrix_gene",
    "plotProfile_gene",
    "plotHeatmap_peak",
    "featureCounts",
    "upset_plot_peaks",
    "chipseeker",
    "maelstrom_report_preparation",
    "infer_strandedness",
    "dupRadar",
    "dupRadar_combine",
    "blind_clustering",
    "merge_volcano_ma",

    "seq2science",
]

rules_to_color = {
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

    # trackhub
    "bam_bigwig": "0.00 0.6 0.85",  # cherry
    "peak_bigpeak": "0.00 0.6 0.85",  # cherry
    "bedgraph_bigwig": "0.00 0.6 0.85",  # cherry
    "trackhub": "0.00 0.6 0.85",  # cherry

    # multiqc
    "multiqc": "0.63 0.6 0.85",  # purple

    # count files
    "coverage_table": "0.28 0.6 0.85",  # green
    "onehot_peaks": "0.28 0.6 0.85",  # green
    "quantile_normalization": "0.28 0.6 0.85",  # green
    "edgeR_normalization": "0.28 0.6 0.85",  # green
    "mean_center": "0.28 0.6 0.85",  # green
    "combine_biological_reps": "0.28 0.6 0.85",  # green

    "gene_id2name": "0.28 0.6 0.85",  # green
    "tpm_matrix": "0.28 0.6 0.85",  # green
    "count_matrix": "0.28 0.6 0.85",  # green
    "txi_count_matrix": "0.28 0.6 0.85",  # green
    "pytxi_count_matrix": "0.28 0.6 0.85",  # green

    "deseq2": "0.28 0.6 0.85",  # green
    "dexseq_count_matrix": "0.28 0.6 0.85",  # green
    "kallistobus_count": "0.28 0.6 0.85",  # green
    "create_bins_SNAP_object": "0.28 0.6 0.85",  # green
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

    def color_node(self, label, color):
        node_id = self.label2id.get(label)
        if node_id is None:
            return

        self.nodes[node_id]["color"] = color

    def remove_node(self, label):
        node_id = self.label2id.get(label)
        if node_id is None:
            return

        # remove the node
        del self.nodes[node_id]
        del self.label2id[label]
        # remove all edges with this node
        for edge in self.edges.copy():
            if node_id in edge:
                self.edges.remove(edge)

    def remove_edge(self, source, target):
        source_id = self.label2id.get(source)
        target_id = self.label2id.get(target)
        edge = (source_id, target_id)
        if edge in self.edges:
            self.edges.remove(edge)

    def hide_node(self, label):
        """remove a node, and connect incoming edges and outgoing edges"""
        node_id = self.label2id.get(label)
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
        self.remove_node(label)
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
