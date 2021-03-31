import re
import os.path
import yaml


def hyperref(text, link):
    if config.get("hyperref", False):
        return f"""<a href="{link}">{text}</a>"""
    else:
        return f"{text} ({link})"


def options(*opts, ret=config):
    for opt in opts:
        assert isinstance(ret, dict), f"Invalid config key '{opt}' in {opts}."
        ret = ret.get(opt)
        if ret is None:
            return ""
    return f" with options '{ret}'"


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
        string = messages[name]

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

        parser = re.compile("(hyperref\(.+?\))")
        while len(parser.findall(string)):
            match = next(parser.finditer(string))
            string = string[:match.span()[0]] + eval(match.group(0)) + string[match.span()[1]:]

        return string


    def text_join(lst, sep=", ", final_sep=" or ", start="", end=""):
        """
        prunes empty strings from a list, joins it with commas and a final separator
        """
        lst = list(filter(None, lst))
        if len(lst) > 1:
            all_but_last = sep.join(lst[:-1])
            last = lst[-1]
            return start + final_sep.join([all_but_last, last]) + end
        return start + "".join(lst) + end

    DESEQ_ANALYSIS = {"rna_seq": "gene expression", "chip_seq": "peak", "atac_seq": "accessibility"}.get(get_workflow(), "")

    messages={
        "bowtie2_align": "Reads were aligned with hyperref('bowtie2 v@bowtie2[bowtie2]','https://dx.doi.org/10.1038%2Fnmeth.1923')" + options("align") + ".",
        "bwa-mem_align": "Reads were aligned with hyperref('bwa-mem v@bwa[bwa]','http://arxiv.org/abs/1303.3997')" + options("align") + ".",
        "bwa-mem2_align": "Reads were aligned with hyperref('bwa-mem2 v@bwamem2[bwa-mem2]','https://arxiv.org/abs/1907.12931')" + options("align") + ".",
        "hisat_splice_aware": "An exon and splice-aware index was generated for HISAT2.",
        "hisat2_align": "Reads were aligned with hyperref('HISAT2 v@hisat2[hisat2]','https://doi.org/10.1038/s41587-019-0201-4')" + options("align") + ".",
        "star_align": "Reads were aligned with hyperref('STAR v@star[star]','https://dx.doi.org/10.1093%2Fbioinformatics%2Fbts635')" + options("align") + ".",
        "sieve_bam":
            text_join(start="Mapped reads were removed if they ",
                      lst=["did not have a minimum mapping quality of {config[min_mapping_quality]}" if config.get("min_mapping_quality", 0) > 0 else "",
                           "were a (secondary) multimapper" if config.get("only_primary_align") else "",
                           "aligned inside the hyperref('ENCODE blacklist','https://doi.org/10.1038/s41598-019-45839-z')" if config.get("remove_blacklist") else "",
                           "had a template length longer than {config[max_template_length]} and shorter than {config[min_template_length]}" if config.get("filter_on_size") else ""],
                      end=" and finally were tn5 bias shifted by seq2science." if config.get("tn5_shift", 0) > 0 else "."),
        # "samtools_sort": "Bam files were sorted with samtools v@samtools[samtools].",
        "sambamba_sort": "Bam files were sorted with hyperref('sambamba v@sambamba[sambamba]','https://doi.org/10.1093/bioinformatics/btv098').",
        "mark_duplicates": "Afterwards, duplicate reads were " + ("removed" if "REMOVE_DUPLICATES=true" in config.get("markduplicates", "") else "marked") + " with hyperref('picard MarkDuplicates v@picard[picard]','http://broadinstitute.github.io/picard').",
        "bam2cram": "Bam files were converted to cram format with samtools v@samtools[samtools].",
        "deseq2":
            text_join(start=f"Differential {DESEQ_ANALYSIS} analysis was performed using hyperref('DESeq2 v@deseq2[bioconductor-deseq2]','https://dx.doi.org/10.1186%2Fs13059-014-0550-8'). To adjust for multiple testing ",
                      lst=[("the (default) Benjamini-Hochberg procedure " if config.get('deseq2', {}).get('multiple_testing_procedure') == "BH" else
                               "hyperref('Independent hypothesis weighting','http://dx.doi.org/10.1038/nmeth.3885') "),
                           "was performed with an FDR cutoff of {config[deseq2][alpha_value]} (default is 0.1). Counts were log transformed using ",
                           "the (default) shrinkage estimator hyperref('apeglm','https://doi.org/10.1093/bioinformatics/bty895'). " if config.get('deseq2', {}).get('shrinkage_estimator') == "apeglm" else (
                               "shrinkage estimator hyperref('ashr','https://doi.org/10.1093/biostatistics/kxw041'). " if config.get('deseq2', {}).get('shrinkage_estimator') == "ashr" else
                                   "the normal prior distribution provided by DESeq2.")], sep=" ", final_sep=" "),
        "count_matrix_txi": "Transcript abundance estimations were aggregated and converted to gene counts using hyperref('tximeta v@tximeta[tximeta]','https://doi.org/10.1101/777888').",
        "run2sra": "Public samples were downloaded from the hyperref('Sequence Read Archive','https://doi.org/10.1093/nar/gkq1019') with help of the ncbi e-utilities.",
        "get_genome": "Genome assembly {wildcards.raw_assembly} was downloaded with hyperref('genomepy {genomepy.__version__}','https://doi.org/10.21105/joss.00320').",
        "custom_extension": "The genome and gene annotations was extended with custom regions.",
        "call_peak_genrich": "Peaks were called with hyperref('genrich v@genrich[genrich]','https://github.com/jsh58/Genrich')" + options("peak_caller", "genrich") + ".",
        "macs2_callpeak": "Peaks were called with hyperref('macs2 v@macs2[macs2]','(https://doi.org/10.1186/gb-2008-9-9-r137')" + options("peak_caller", "macs2") + " in {params.format} mode. The effective genome size was estimated by taking the number of unique kmers in the assembly of the same length as the average read length for each sample.",
        "keep_mates": "After alignment we removed paired-end info from reads with seq2science to utilize both mates in the paired-end reads.",
        "idr": "Narrowpeak files of biological replicates belonging to the same condition were merged with the hyperref('irreproducible discovery rate v@idr[idr]','http://dx.doi.org/10.1214/11-AOAS466').",
        "macs_cmbreps": "Narrowpeak files of biological replicates belonging to the same condition were merged with fisher's method in macs2 v@macs2[macs2].",
        "multiqc": "Quality control metrics were aggregated by hyperref('MultiQC v@multiqc[multiqc]','http://dx.doi.org/10.1093/bioinformatics/btw354').",
        "samtools_stats": "General alignment statistics were collected by hyperref('samtools stats v@samtools[samtools]','https://doi.org/10.1093/bioinformatics/btp352').",
        "featureCounts_qc": "The fraction reads in peak score (frips) was calculated by hyperref('featurecounts v@subread[subread]','https://doi.org/10.1093/bioinformatics/btt656').",
        "fastqc": "Fastq quality was measured by hyperref('FastQC v@fastqc[fastqc]','http://www.bioinformatics.babraham.ac.uk/projects/fastqc').",
        "computeMatrix": "hyperref('Deeptools v@deeptools[deeptools]','https://doi.org/10.1093/nar/gkw257') was used for the fingerprint, profile, correlation and dendrogram/heatmap plots, where the heatmap was made" + options("deeptools_multibamsummary") + ".",
        "dupradar": "RNA-seq read duplication types were analyzed using hyperref('dupRadar v@dupradar[bioconductor-dupradar]','https://doi.org/10.1186/s12859-016-1276-2').",
        "decoy_transcripts": "Decoy transcript were generated in order improve improve Salmon indexing accuracy (using the script from https://github.com/COMBINE-lab/SalmonTools)",
        "salmon_quant": "Transcript abundances were quantified with hyperref('Salmon v@salmon[salmon]','https://doi.org/10.1038/nmeth.4197')" + options("quantifier_flags") + ".",
        "htseq_count": "Read counting and summarizing to gene-level was performed on filtered bam using hyperref('HTSeq-count v@gene_counts[htseq]','https://doi.org/10.1093/bioinformatics/btu638').",
        "kallistobus-count": "Reads were aligned and transformed to bus format with kb-python v@kallistobus[kb-python], a python wrapper for hyperref('kallisto','https://doi.org/10.1038/nbt.3519') and hyperref('bustools','doi.org/10.1101/673285').",
        "featurecounts_rna": "Read counting and summarizing to gene-level was performed on filtered bam using hyperref('featureCounts v@gene_counts[subread]','https://doi.org/10.1093/bioinformatics/btt656').",
        "dexseq": "Additionally, exon usage was counted using hyperref('[DEXSeq] v@dexseq[bioconductor-dexseq]','https://doi.org/doi:10.18129/B9.bioc.DEXSeq') for (potential) downstream analysis.",
        "create_bins_SNAP_object": "We used hyperref('snaptools v@snaptools[snaptools]','https://doi.org/10.1101/615179') to create a snapobject" + options("snaptools_opt") + " and added a binned genome matrix" + options("bin_opt") + ".",
        "infer_strandedness": "Sample sequencing strandedness was inferred using hyperref('RSeQC v@gene_counts[rseqc]','https://doi.org/10.1093/bioinformatics/bts356') in order to improve quantification accuracy.",
        "trackhub": "We used the hyperref('UCSC genome browser','http://www.genome.org/cgi/doi/10.1101/gr.229102') to visualize and inspect alignment.",
        "trimgalore_SE": "We trimmed single-end reads with hyperref('trim galore! v@trimgalore[trim-galore]','http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/')" + options("trimoptions") + " and hyperref('cutadapt','https://doi.org/10.14806/ej.17.1.200').",
        "trimgalore_PE": "We trimmed paired-end reads with hyperref('trim galore! v@trimgalore[trim-galore]','http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/')" + options("trimoptions") + " and hyperref('cutadapt','https://doi.org/10.14806/ej.17.1.200').",
        "fastp_SE": "We trimmed single-end reads with hyperref('fastp v@fastp[fastp]','https://doi.org/10.1093/bioinformatics/bty560')" + options("trimoptions") + ".",
        "fastp_PE": "We trimmed paired-end reads with hyperref('fastp v@fastp[fastp]','https://doi.org/10.1093/bioinformatics/bty560')" + options("trimoptions") + ".",
        "chipseeker": "A peak feature distribution plot and peak localization plot relative to TSS were made with hyperref('chipseeker','https://doi.org/doi:10.18129/B9.bioc.ChIPseeker').",  # v@chipseeker[chipseeker]
        "combine_peaks": "A consensus set of summits was made with hyperref('gimmemotifs.combine_peaks v@gimme[gimmemotifs]','https://www.biorxiv.org/content/10.1101/474403v1.full').",
        "bed_slop": "All summits were extended with '{config[slop]}' to get a consensus peakset.",
        "coverage_table": "And finally we made a count table from the conensus peakset with hyperref('gimmemotifs.combine_peaks v@gimme[gimmemotifs]','https://www.biorxiv.org/content/10.1101/474403v1.full').",
    }
