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
scripts_dir    <- snakemake@params$scripts_dir
assembly       <- snakemake@wildcards$assembly
out_png        <- snakemake@output[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# load custom functions
source(file.path(scripts_dir, "utils.R"))

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('threads      <- ',  threads,      '\n')
cat('counts_file  <- "', counts_file,  '"\n', sep = "")
cat('samples_file <- "', samples_file, '"\n', sep = "")
cat('replicates   <- ',  replicates,   '\n')
cat('assembly     <- "', assembly,     '"\n', sep = "")
cat('out_png      <- "', out_png,      '"\n', sep = "")
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')


## obtain coldata, the metadata input for DESeq2
samples <- parse_samples(samples_file, assembly, replicates)
# samples <- read.delim(samples_file, sep = "\t", na.strings = "", comment.char = "#", stringsAsFactors = F)
# if ("technical_replicate" %in% colnames(samples) & isTRUE(replicates)) {
#   samples$technical_replicate[is.na(samples$technical_replicate)] <- as.character(samples$sample[is.na(samples$technical_replicate)])
#   samples <- subset(samples, !duplicated(technical_replicate))
#   row.names(samples) <- samples$technical_replicate
# } else {
#   row.names(samples) <- samples$sample
# }

# filter for assembly, remove NAs and add random variables (not needed for blind clustering)
cols <- ifelse("descriptive_name" %in% colnames(samples), c('assembly', 'descriptive_name'), 'assembly')
coldata  <- samples[cols]
coldata['assembly'] <- factor(as.character(seq_len(nrow(coldata))))
# # filter for assembly, remove NAs and add random variables (not needed for blind clustering)
# cols <- ifelse("descriptive_name" %in% colnames(samples), list('assembly', 'descriptive_name'), list('assembly'))
# coldata  <- samples[samples$assembly == assembly, cols, drop = F]

# filter counts to speed up DESeq
counts <- read.table(counts_file, row.names = 1, header = T, stringsAsFactors = F, sep = '\t', check.names = F)
reduced_counts <- counts[rowSums(counts) > 0, colnames(counts) %in% rownames(coldata)]


## DESeq2
dds <- run_deseq2(reduced_counts, coldata, ~ 1, threads)
# # setup parallelization
# parallel <- FALSE
# if (threads > 1) {
#   register(MulticoreParam(threads))
#   parallel <- TRUE
# }
#
#
# # Combine counts and metadata
# dds <- DESeqDataSetFromMatrix(countData = reduced_counts,
#                               colData = coldata,
#                               design = ~ 1) #dummy design, we only want the sample correlations
#
# # normalization and preprocessing
# dds <- DESeq(dds, parallel=parallel)
# cat('\n')

## Heatmaps
# list of aesthetics for all heatmaps
heatmap_aes <- heatmap_aesthetics(nrow(coldata))

# transform counts for all heatmaps
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsdMatrix <- t(assay(vsd))

# Heatmap of the sample-to-sample distances
# see http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-sample-to-sample-distances
mat <- as.matrix(dist(vsdMatrix))
heatmap_names(distMatrix, coldata)

title <- 'Sample distance clustering (blind)'
legend_aes <- list(
  "breaks" = c(min(sampleDistMatrix), mean(c(min(sampleDistMatrix), max(sampleDistMatrix))), max(sampleDistMatrix)),
  "labels" = c("low",                 "distance\n(euclidean)",                               "high"),
)
out_pdf <- sub(".png", ".pdf", out_png)
heatmap_plot(mat, title, heatmap_aes, legend_aes, out_pdf)

# convert pdf to png (required for MULTIQC)
pdf2png(out_pdf, out_png)


# Heatmaps of various correlation matrices
for (method in list("spearman", "pearson")){
  # Kendall is also available, but is computationally much more complex
  mat <- as.matrix(cor(vsdMatrix, method=method))
  heatmap_names(mat, coldata)

  # capitalize method name
  name <- paste0(toupper(substring(method, 1, 1)), substring(method, 2))
  title <- paste(name, 'correlation clustering')
  legend_aes <- list(
    "breaks" = NA,
    "labels" = NA,
  )
  out_pdf_cor <- sub("Sample_distance_clustering", gsub(" ", "_", title), out_pdf)
  heatmap_plot(mat, title, heatmap_aes, legend_aes, out_pdf_cor)

  # convert pdf to png (required for MULTIQC)
  out_png_cor <- sub(".pdf", ".png", out_pdf_cor)
  pdf2png(out_pdf_cor, out_png_cor)
}


# spearmanMatrix <- as.matrix(cor(vsdMatrix, method="spearman"))
# pearsonMatrix  <- as.matrix(cor(vsdMatrix, method="pearson"))
# kendallMatrix  <- as.matrix(cor(vsdMatrix, method="kendall"))



# # transform the data
# main <- 'Sample distance clustering (blind)'
# log_counts <- varianceStabilizingTransformation(dds, blind = TRUE)
# sampleDistMatrix <- as.matrix(dist(t(assay(log_counts))))
# rownames(sampleDistMatrix) <- if (length(cols) == 1) {colnames(counts(dds))} else {coldata$descriptive_name}
# colnames(sampleDistMatrix) <- if (length(cols) == 1) {colnames(counts(dds))} else {coldata$descriptive_name}
#
#
#
#
#
# # if (num_samples < 16) {
# #   cell_dimensions <- as.integer(160/num_samples)  # minimal size to fit the legend
# #   fontsize        <- 8
# # } else if (num_samples < 24) {
# #   cell_dimensions <- 10  # pleasant size
# #   fontsize        <- 8
# # } else if (num_samples < 32) {
# #   cell_dimensions <- as.integer(25 - 0.625*num_samples)  # linear shrink
# #   fontsize        <- 8 - 0.15*num_samples
# # } else {
# #   cell_dimensions <- 5  # minimal size
# #   fontsize        <- 3.2
# # }
# # show_colnames <- ifelse(num_samples > 28, TRUE, FALSE)
# # show_rownames <- ifelse(num_samples > 28, FALSE, TRUE)
#
# #
# # cell_dimensions <- (if (num_samples < 16) {as.integer(160/num_samples)}             # minimal size to fit the legend
# #                    else if (num_samples < 24) {10}                                  # pleasant size
# #                    else if (num_samples < 32) {as.integer(25 - 0.625*num_samples)}  # linear shrink
# #                    else {5})                                                        # minimal size
# # fontsize <- (if (num_samples < 16) {8}
# #             else if (num_samples < 32) {8 - 0.15*num_samples}
# #             else {3.2})
#
# # make heatmap and save as pdf (only pdfs can consistently save text properly)
# out_pdf <- sub(".png", ".pdf", out_png)
# pheatmap(
#   sampleDistMatrix,
#   main = main,
#   angle_col = 45,
#   show_colnames = if (num_samples > 28) {TRUE} else {FALSE},  # show names underneath if the image gets to wide
#   show_rownames = if (num_samples > 28) {FALSE} else {TRUE},
#   fontsize = fontsize,
#   legend_breaks = c(min(sampleDistMatrix), mean(c(min(sampleDistMatrix), max(sampleDistMatrix))), max(sampleDistMatrix)),
#   legend_labels = c("low",                 "distance\n(euclidean)",                               "high"),
#   cellwidth  = cell_dimensions,
#   cellheight = cell_dimensions,
#   col=colors,
#   filename=out_pdf,
# )
#
# # convert pdf to png (required for MULTIQC)
# pdf2png(out_pdf, out_png)
# # pdftools::pdf_convert(
# #   out_pdf,
# #   format = "png",
# #   filenames = out_png,
# #   dpi = 300,  # standard minimum
# #   antialias = TRUE,
# #   verbose = TRUE
# # )
