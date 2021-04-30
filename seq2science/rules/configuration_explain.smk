import os.path
import yaml

if not config.get("explain_rule"):
    def explain_rule(name):
        """
        Does not return a string explaining the workflow by parsing all messages
        """
        return None

else:
    REFERENCES = {
        'apeglm': 'https://doi.org/10.1093/bioinformatics/bty895',
        'ashr': 'https://doi.org/10.1093/biostatistics/kxw041',
        'bowtie2': 'https://dx.doi.org/10.1038%2Fnmeth.1923',
        'bustools': 'https://doi.org/10.1101/673285',
        'bwa': 'http://arxiv.org/abs/1303.3997',
        'bwa-mem2': 'https://arxiv.org/abs/1907.12931',
        'chipseeker': 'https://doi.org/doi:10.18129/B9.bioc.ChIPseeker',
        'cutadapt': 'https://doi.org/10.14806/ej.17.1.200',
        'deeptools': 'https://doi.org/10.1093/nar/gkw257',
        'deseq2': 'https://dx.doi.org/10.1186%2Fs13059-014-0550-8',
        'dexseq': 'https://doi.org/doi:10.18129/B9.bioc.DEXSeq',
        'dupradar': 'https://doi.org/10.1186/s12859-016-1276-2',
        'encode blacklist': 'https://doi.org/10.1038/s41598-019-45839-z',
        'fastp': 'https://doi.org/10.1093/bioinformatics/bty560',
        'fastqc': 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc',
        'genomepy': 'https://doi.org/10.21105/joss.00320',
        'genrich': 'https://github.com/jsh58/Genrich',
        'gimmemotifs': 'https://www.biorxiv.org/content/10.1101/474403v1.full',
        'hisat2': 'https://doi.org/10.1038/s41587-019-0201-4',
        'htseq': 'https://doi.org/10.1093/bioinformatics/btu638',
        'idr': 'http://dx.doi.org/10.1214/11-AOAS466',
        'independent hypothesis weighting': 'http://dx.doi.org/10.1038/nmeth.3885',
        'kallisto': 'https://doi.org/10.1038/nbt.3519',
        'macs2': 'https://doi.org/10.1186/gb-2008-9-9-r137',
        'multiqc': 'http://dx.doi.org/10.1093/bioinformatics/btw354',
        'picard': 'http://broadinstitute.github.io/picard',
        'rseqc': 'https://doi.org/10.1093/bioinformatics/bts356',
        'salmon': 'https://doi.org/10.1038/nmeth.4197',
        'sambamba': 'https://doi.org/10.1093/bioinformatics/btv098',
        'samtools': 'https://doi.org/10.1093/bioinformatics/btp352',
        'sequence read archive': 'https://doi.org/10.1093/nar/gkq1019',
        'snaptools': 'https://doi.org/10.1101/615179',
        'star': 'https://dx.doi.org/10.1093%2Fbioinformatics%2Fbts635',
        'subread': 'https://doi.org/10.1093/bioinformatics/btt656',
        'trim-galore': 'http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/',
        'tximeta': 'https://doi.org/10.1101/777888',
        'ucsc genome browser': 'http://www.genome.org/cgi/doi/10.1101/gr.229102',
        'seurat': 'https://satijalab.org/seurat/',
        'kb_seurat_pp': 'https://github.com/Rebecza/scRNA-seq',
    }

    def explain_rule(name):
        """
        Return a string explaining the workflow by parsing all messages
        """
        return MESSAGES[name]


    ENV_DIR = os.path.normpath(os.path.join(config['rule_dir'],"..","envs"))

    def version(tool, env=None):
        """
        Return the tool version from the env.
        Returns an empty string if no version can be found.
        """
        tool = tool.lower()  # conda always uses lower case
        env = env if env else tool  # most tools are in a similarly named env
        yaml_file = f"{ENV_DIR}/{env}.yaml"
        with open(yaml_file,'r') as stream:
            env = yaml.safe_load(stream)

        for dependency in env["dependencies"]:
            if tool in dependency:
                return f" v{dependency[dependency.find('=') + 1:]}"
        return ""


    def hyperref(tool, text=None):
        text = text if text else tool
        ref = REFERENCES[tool.lower()]
        if config.get("hyperref",False):
            return f"""<a href="{ref}">{text}</a>"""
        return f"{text} ({ref})"


    def href_v(tool, text=None, env=None):
        text = text if text else tool
        v = version(tool,env)
        return hyperref(tool,text + v)


    def text_join(lst, sep=", ", final_sep=" or ", start="", end=""):
        """
        prunes empty strings from a list, joins it with commas and a final separator
        """
        lst = list(filter(None,lst))
        if len(lst) > 1:
            all_but_last = sep.join(lst[:-1])
            last = lst[-1]
            return start + final_sep.join([all_but_last, last]) + end
        return start + "".join(lst) + end


    def options(*opts, ret=config):
        for opt in opts:
            assert isinstance(ret,dict), f"Invalid config key '{opt}' in {opts}."
            ret = ret.get(opt)
            if ret in [None, ""]:
                return ""
        return f" with options '{ret}'"


    DESEQ_ANALYSIS = {"rna_seq": "gene expression", "chip_seq": "peak", "atac_seq": "accessibility"}.get(get_workflow(),"")

    MESSAGES = {
        "bowtie2_align": f"Reads were aligned with {href_v('bowtie2')}{options('align')}.",
        "bwa-mem_align": f"Reads were aligned with {href_v('bwa',text='bwa-mem')}{options('align')}.",
        "bwa-mem2_align": f"Reads were aligned with {href_v('bwa-mem2',env='bwamem2')}{options('align')}.",
        "hisat_splice_aware": "An exon and splice-aware index was generated for HISAT2.",
        "hisat2_align": f"Reads were aligned with {href_v('HISAT2')}{options('align')}.",
        "star_align": f"Reads were aligned with {href_v('STAR')}{options('align')}.",
        "sieve_bam":
            text_join(start="Mapped reads were removed if they ",
                lst=[f"did not have a minimum mapping quality of {config['min_mapping_quality']}" if config.get("min_mapping_quality",0) > 0 else "",
                     "were a (secondary) multimapper" if config.get("only_primary_align") else "",
                     f"aligned inside the {hyperref('ENCODE blacklist')}" if config.get("remove_blacklist") else "",
                     f"had a template length longer than {config['max_template_length']} bp and shorter than {config['min_template_length']} bp" if config.get("filter_on_size") else ""],
                end=" and finally were tn5 bias shifted by seq2science." if config.get("tn5_shift",0) > 0 else "."),
        "sambamba_sort": f"Bam files were sorted with {href_v('sambamba')}.",
        "mark_duplicates": f"Afterwards, duplicate reads were {'removed' if 'REMOVE_DUPLICATES=true' in options('markduplicates') else 'marked'} with {href_v('picard',text='Picard MarkDuplicates')}.",
        "bam2cram": f"Bam files were converted to cram format with samtools {href_v('samtools')}.",
        "deseq2": "" if (deseq_opts := config.get('deseq2')) is None else \
            text_join(
                start=f"Differential {DESEQ_ANALYSIS} analysis was performed using {href_v('DESeq2')}. To adjust for multiple testing ",
                lst=["the (default) Benjamini-Hochberg procedure " if deseq_opts['multiple_testing_procedure'] == "BH" else href_v('Independent hypothesis weighting',env='deseq2'),
                     f"was performed with an FDR cutoff of {deseq_opts['alpha_value']} (default is 0.1). Counts were log transformed using ",
                     f"the (default) shrinkage estimator {href_v('apeglm',env='deseq2')}. " if deseq_opts['shrinkage_estimator'] == "apeglm" else (
                         f"shrinkage estimator {href_v('ashr',env='deseq2')}. " if deseq_opts['shrinkage_estimator'] == "ashr" else
                            "the normal prior distribution provided by DESeq2.")],sep=" ",final_sep=" "),
        "count_matrix_txi": f"Transcript abundance estimations were aggregated and converted to gene counts using {href_v('tximeta')}.",
        "run2sra": f"Public samples were downloaded from the {hyperref('Sequence Read Archive')} with help of the ncbi e-utilities.",
        "get_genome": f"Genome assembly {{wildcards.raw_assembly}} was downloaded with {hyperref('genomepy',text=f'genomepy {genomepy.__version__}')}.",
        "custom_extension": "The genome and gene annotations was extended with custom regions.",
        "call_peak_genrich": f"Peaks were called with {href_v('genrich')}{options('peak_caller','genrich')}.",
        "macs2_callpeak": f"Peaks were called with {href_v('macs2')}{options('peak_caller','macs2')} in {{params.format}} mode. The effective genome size was estimated by taking the number of unique kmers in the assembly of the same length as the average read length for each sample.",
        "keep_mates": "After alignment we removed paired-end info from reads with seq2science to utilize both mates in the paired-end reads.",
        "idr": f"Narrowpeak files of biological replicates belonging to the same condition were merged with the {href_v('idr',text='irreproducible discovery rate')}.",
        "macs_cmbreps": "Narrowpeak files of biological replicates belonging to the same condition were merged with fisher's method in macs2.",
        "multiqc": f"Quality control metrics were aggregated by {href_v('MultiQC')}.",
        "samtools_stats": f"General alignment statistics were collected by {href_v('samtools',text='samtools stats')}.",
        "featureCounts_qc": f"The fraction reads in peak score (frips) was calculated by {href_v('subread',text='featurecounts')}.",
        "fastqc": f"Fastq quality was measured by {href_v('FastQC')}.",
        "computeMatrix": f"{href_v('Deeptools')} was used for the fingerprint, profile, correlation and dendrogram/heatmap plots{', where the heatmap was made' + dt_opts if (dt_opts := options('deeptools_multibamsummary')) else ''}.",
        "dupradar": f"RNA-seq read duplication types were analyzed using {href_v('dupRadar')}.",
        "decoy_transcripts": "Decoy transcript were generated in order improve improve Salmon indexing accuracy using the scripts provided in the Salmon manual",
        "salmon_quant": f"Transcript abundances were quantified with {href_v('Salmon')}{options('quantifier_flags')}.",
        "htseq_count": f"Read counting and summarizing to gene-level was performed on filtered bam using {href_v('htseq',text='HTSeq-count',env='gene_counts')}.",
        "kallistobus-count": f"Reads were aligned and transformed to bus format with kb-python{version('kb-python','kallistobus')}, a python wrapper for {hyperref('kallisto')} and {hyperref('bustools')}.",
        "featurecounts_rna": f"Read counting and summarizing to gene-level was performed on filtered bam using {href_v('subread',text='featureCounts',env='gene_counts')}.",
        "dexseq": f"Additionally, exon usage was counted using {href_v('DEXSeq')} for (potential) downstream analysis.",
        "create_bins_SNAP_object": f"We used {href_v('snaptools')} to create a snapobject{options('snaptools_opt')} and added a binned genome matrix{options('bin_opt')}.",
        "infer_strandedness": f"Sample sequencing strandedness was inferred using {href_v('RSeQC',env='gene_counts')} in order to improve quantification accuracy.",
        "trackhub": f"We used the {hyperref('UCSC genome browser')} to visualize and inspect alignment.",
        "trimgalore_SE": f"We trimmed single-end reads with {href_v('trim-galore',text='trim galore!',env='trimgalore')}{options('trimoptions')} and {hyperref('cutadapt')}.",
        "trimgalore_PE": f"We trimmed paired-end reads with {href_v('trim-galore',text='trim galore!',env='trimgalore')}{options('trimoptions')} and {hyperref('cutadapt')}.",
        "fastp_SE": f"We trimmed single-end reads with {href_v('fastp')}{options('trimoptions')}.",
        "fastp_PE": f"We trimmed paired-end reads with {href_v('fastp')}{options('trimoptions')}.",
        "chipseeker": f"A peak feature distribution plot and peak localization plot relative to TSS were made with {hyperref('chipseeker')}.",  # TODO: replace with href_v
        "combine_peaks": f"A consensus set of summits was made with {href_v('gimmemotifs',text='gimmemotifs.combine_peaks',env='gimme')}.",
        "bed_slop": f"All summits were extended with {config.get('slop')} bp to get a consensus peakset.",
        "coverage_table": f"And finally we made a count table from the conensus peakset with gimmemotifs.",  # already cited in "combine_peaks"
        "kb_seurat_pp": f"scRNA count post-processing was performed using the {hyperref('seurat') } based {hyperref('kb_seurat_pp')} markdown with the following settings{options('seurat')}.",
    }
