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
replicates     <- snakemake@params$replicates
assembly       <- snakemake@wildcards$assembly
out_plot       <- snakemake@output[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('threads      <-',   threads, '\n')
cat('log_file     <- "', log_file,     '"\n', sep = "")
cat('counts_file  <- "', counts_file,  '"\n', sep = "")
cat('samples_file <- "', samples_file, '"\n', sep = "")
cat('replicates   <- ',  replicates, '\n')
cat('assembly     <- "', assembly,     '"\n', sep = "")
cat('out_plot     <- "', out_plot,     '"\n', sep = "")
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')


## obtain coldata, the metadata input for DESeq2
samples <- read.delim(samples_file, sep = "\t", na.strings = "", comment.char = "#", stringsAsFactors = F)
if ("replicate" %in% colnames(samples) & isTRUE(replicates)) {
  samples$replicate[is.na(samples$replicate)] <- as.character(samples$sample[is.na(samples$replicate)])
  samples <- subset(samples, !duplicated(replicate))
  row.names(samples) <- samples$replicate
} else {
  row.names(samples) <- samples$sample
}

# filter for assembly, remove NAs and add random variables (not needed for blind clustering)
coldata  <- samples[samples$assembly == assembly, 'assembly', drop = F]
coldata[,1] <- factor(as.character(c(1:nrow(coldata))))


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
log_counts <- varianceStabilizingTransformation(dds, blind = TRUE)
sampleDistMatrix <- as.matrix(dist(t(assay(log_counts))))
rownames(sampleDistMatrix) <- colnames(counts(dds))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

svg(file=out_plot)
pheatmap(sampleDistMatrix,
         main = main,
         angle_col = 45,
         col=colors)
invisible(dev.off())
