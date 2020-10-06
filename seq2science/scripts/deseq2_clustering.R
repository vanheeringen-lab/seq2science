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
cols = if ("descriptive_name" %in% colnames(samples)) {c('assembly', 'descriptive_name')} else {c('assembly')}
coldata  <- samples[samples$assembly == assembly, cols, drop = F]
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
rownames(sampleDistMatrix) <- if (length(cols) == 1) {colnames(counts(dds))} else {coldata$descriptive_name}
colnames(sampleDistMatrix) <- if (length(cols) == 1) {colnames(counts(dds))} else {coldata$descriptive_name}
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

num_samples <- dim(coldata)[1]
cell_dimensions = (if (num_samples < 16) {as.integer(160/num_samples)}              # minimal size to fit the legend
                   else if (num_samples < 24) {10}                                  # pleasant size
                   else if (num_samples < 32) {as.integer(25 - 0.625*num_samples)}  # linear shrink
                   else {5})                                                        # minimal size
pheatmap(sampleDistMatrix,
         main = main,
         angle_col = 45,
         show_colnames = if (num_samples > 28) {TRUE} else {FALSE},  # show names underneath if the image gets to wide
         show_rownames = if (num_samples > 28) {FALSE} else {TRUE},
         fontsize = 8,
         legend_breaks = c(min(sampleDistMatrix), max(sampleDistMatrix)),
         legend_labels = c("high", "low"),
         cellwidth  = cell_dimensions,
         cellheight = cell_dimensions,
         col=colors,
         filename=out_plot
)
