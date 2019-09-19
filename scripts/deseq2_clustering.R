suppressMessages({
  library(DESeq2)
  library(BiocParallel)
  library(RColorBrewer)
  library(pheatmap)
})

# snakemake variables
threads        <- snakemake@threads[[1]]
log_file       <- snakemake@log[[1]]
counts_file    <- snakemake@input[[1]]
samples_file   <- snakemake@params$samples
contrast       <- snakemake@params$some_contrast
assembly       <- snakemake@wildcards$assembly
out_plot       <- snakemake@output[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")


## parse the design contrast
# a contrast is always in the form 'batch+condition_group1_group2', where batch(+) is optional
# extracting batch, condition and groups
contr <- gsub('~' ,'' , contrast)
if (grepl('\\+', contr)) {
  contr <- strsplit(contr, '\\+')[[1]][2]
}
condition <- strsplit(contr, '_')[[1]][1]


## obtain coldata, the metadata input for DESeq2
samples <- read.delim(samples_file, row.names=1, na.strings = "")
# filter for assembly and remove NAs
coldata  <- samples[samples$assembly == assembly, condition, drop = F]


## filter counts to speed up DESeq
counts <- read.table(counts_file, row.names = 1, header = T, stringsAsFactors = F, sep = '\t', check.names = F)
reduced_counts <- counts[rowSums(counts) > 0, colnames(counts) %in% rownames(coldata)]


## DESeq2
# setup parallelization
parallel <- FALSE
if (threads > 1) {
  register(MulticoreParam(threads))
  parallel <- TRUE
}


# Combine counts and metadata
dds <- DESeqDataSetFromMatrix(countData = reduced_counts,
                              colData = coldata,
                              design = ~ 1) #dummy design, we only want the sample correlations

# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)
cat('\n')

# transform the data
main = 'Sample distance clustering (blind)'
log_counts <- vst(dds, blind = TRUE)
sampleDistMatrix <- dist(t(assay(log_counts)))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

svg(file=out_plot)
pheatmap(sampleDistMatrix,
         main = main,
         col=colors)
invisible(dev.off())
