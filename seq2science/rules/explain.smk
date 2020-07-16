import re
import os.path
import yaml


if not config.get("explain_rule"):
    def explain_rule(name):
        """
        Does not return a string explaining the workflow by parsing all messages
        """
        return None

else:
    def explain_rule(name):
        """
        Return a string explaining the workflow by parsing all messages
        """
        try:
            string = messages[name]
        except:
            string = name

        # clean our explanation
        string = string.replace("\n", "")
        string = " ".join(string.split())

        # find our environment
        env_dir = os.path.normpath(os.path.join(config['rule_dir'], "..", "envs"))

        parser = re.compile("@([^[]*)\[([^]]*)\]")
        while len(parser.findall(string)):
            match = next(parser.finditer(string))

            # parse our environment
            yaml_file = f"{env_dir}/{match.group(1)}.yaml"
            tool = match.group(2)
            with open(yaml_file, 'r') as stream:
                env = yaml.safe_load(stream)

            for dependency in env["dependencies"]:
                if tool in dependency:
                    version = dependency[dependency.find("=") + 1:]
                    break
            else:
                continue

            # replace the placeholder with the actual version
            string = string[:match.span()[0]] + version + string[match.span()[1]:]
        return string

    def text_join(lst, final_sep=" or ", start="", end=""):
        """
        prunes empty strings from a list, joins it with commas and a final separator
        """
        lst = list(filter(None, lst))
        if len(lst) > 1:
            all_but_last = ', '.join(lst[:-1])
            last = lst[-1]
            return start + final_sep.join([all_but_last, last]) + end
        return start + lst[0] + end

    messages={
        "bowtie2_align":"Reads were aligned with bowtie2 v@bowtie2[bowtie2] (https://dx.doi.org/10.1038%2Fnmeth.1923) with options '{config[align]}'.",
        "bwa_mem": "Reads were aligned with bwa-mem v@bwa[bwa] (http://arxiv.org/abs/1303.3997) with options '{config[align]}'.",
        "bwa_mem2": "TODO",
        "hisat2_align": "Reads were aligned with hisat2 v@hisat2[hisat2] (https://doi.org/10.1038/s41587-019-0201-4) with options '{config[align]}'.",
        "star_align": "Reads were aligned with STAR v@star[star] () with options '{config[align]}'.",
        "sieve_bam":
            text_join(start="Mapped reads were removed if they ",
                      lst=["did not have a minimum mapping quality of {config[min_mapping_quality]}" if config.get("min_mapping_quality", 0) > 0 else "",
                           "were a (secondary) multimapper" if config.get("only_primary_align") else "",
                           "aligned inside the ENCODE blacklist (https://doi.org/10.1038/s41598-019-45839-z)" if config.get("remove_blacklist") else ""],
                      end=" and finally were tn5 bias shifted by seq2science." if config.get("tn5_shift", 0) > 0 else "."),
        "samtools_sort": "Bam files were sorted with samtools v@samtools[samtools].",
        "sambamba_sort": "Bam files were sorted with sambamba v@sambamba[sambamba].",
        "mark_duplicates": "Afterwards duplicate reads were marked with picard MarkDuplicates v@picard[picard].",
        "bam2cram": "Bam files were converted to cram format with samtools v@samtools[samtools].",
        "deseq2": "TODO",
        "blind_clustering": "TODO",
        "count_matrix_txi": "Transcript abundance estimations were aggregated and converted to gene counts using tximeta v@tximeta[tximeta]",
        "id2sra": "Public samples were downloaded from the Sequence Read Archive (https://doi.org/10.1093/nar/gkq1019) with help of the ncbi e-utilities.",
        "get_genome": "Genomic assembly {wildcards.assembly} was downloaded with genomepy v@get_genome[genomepy] (https://doi.org/10.21105/joss.00320).",
        "call_peak_genrich":"Peaks were called with genrich v@genrich[genrich] (https://github.com/jsh58/Genrich) with options '{config[peak_caller][genrich]}'.",
        "macs2_callpeak": "Peaks were called with macs2 v@macs2[macs2] with options '{config[peak_caller][macs2]}' in {params.format} mode.",
        "keep_mates": "After alignment we removed paired-end info from reads with seq2science to utilize both mates in the paired-end reads.",
        "idr": "Narrowpeak files of replicates belonging to the same condition were merged with the irreproducible discovery rate v@idr[idr].",
        "macs_cmbreps": "Narrowpeak files of replicates belonging to the same condition were merged with fisher's method in macs2 v@macs2[macs2].",
        "multiqc": "Quality control metrics were aggregated by MultiQC v@qc[multiqc].",
        "samtools_stats": "General alignment statistics were collected by samtools stats v@samtools[samtools].",
        "featureCounts_qc":"The fraction reads in peak score (frips) was calculated by featurecounts v@subread[subread] (https://doi.org/10.1093/bioinformatics/btt656).",
        "fastqc": "Fastq quality was measured by FastQC v@qc[fastqc].",
        "computeMatrix": "Deeptools v@deeptools[deeptools] (https://doi.org/10.1093/nar/gkw257) was used for the fingerprint, profile, correlation and heatmap plots.",
        "decoy_transcripts": "Decoy transcript were generated in order improve improve Salmon indexing accuracy (using the script from https://github.com/COMBINE-lab/SalmonTools)",
        "salmon_quant": "TODO",
        "htseq_count":"TODO",
        "featurecounts_rna":"TODO",
        "create_bins_SNAP_object": "We used snaptools v@snaptools[snaptools] to create a snapobject with options config[snaptools_opt] and added a binned genome matrix with options {config[bin_opt]}.",
        "trackhub": "We used the UCSC genome browser (http://www.genome.org/cgi/doi/10.1101/gr.229102) to visualize and inspect alignment.",
        "trim_galore_SE": "We trimmed single-end reads with trim galore! v@trimgalore[trim-galore] (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) with options '{config[trim_galore]}' and cutadapt (https://doi.org/10.14806/ej.17.1.200).",
        "trim_galore_PE": "We trimmed paired-end reads with trim galore! v@trimgalore[trim-galore] (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) with options '{config[trim_galore]}' and cutadapt (https://doi.org/10.14806/ej.17.1.200).",
    }



    # def explain_rule(string):
    #     """
    #     Parse a message.
    #     """
    #     # clean our explanation
    #     string = string.replace("\n", "")
    #     string = " ".join(string.split())
    #
    #     # find our environment
    #     env_dir = os.path.normpath(os.path.join(config['rule_dir'], "..", "envs"))
    #
    #     parser = re.compile("@([^[]*)\[([^]]*)\]")
    #     while len(parser.findall(string)):
    #         match = next(parser.finditer(string))
    #
    #         # parse our environment
    #         yaml_file = f"{env_dir}/{match.group(1)}.yaml"
    #         tool = match.group(2)
    #         with open(yaml_file, 'r') as stream:
    #             env = yaml.safe_load(stream)
    #
    #         for dependency in env["dependencies"]:
    #             if tool in dependency:
    #                 version = dependency[dependency.find("=") + 1:]
    #                 break
    #         else:
    #             continue
    #
    #         # replace the placeholder with the actual version
    #         string = string[:match.span()[0]] + version + string[match.span()[1]:]
    #     return string
