suppressMessages({
  library(DESeq2)
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

# filter for assembly, remove NAs and add random variables (not needed for blind clustering)
cols <- if("descriptive_name" %in% colnames(samples)){c('assembly', 'descriptive_name')}else{'assembly'}
coldata  <- samples[cols]
coldata['assembly'] <- factor(as.character(seq_len(nrow(coldata))))

# filter counts to speed up DESeq
counts <- read.table(counts_file, row.names = 1, header = T, stringsAsFactors = F, sep = '\t', check.names = F)
# has_descriptive <- "descriptive_name" %in% colnames(coldata)
# names <- if (has_descriptive) {coldata$descriptive_name} else {rownames(coldata)}
# reduced_counts <- counts[rowSums(counts) > 0, colnames(counts) %in% names]
reduced_counts <- counts[rowSums(counts) > 0, rownames(coldata)]


## DESeq2
dds <- run_deseq2(reduced_counts, coldata, ~ 1, threads)


## Heatmaps
# list of aesthetics for all heatmaps
heatmap_aes <- heatmap_aesthetics(nrow(coldata))

# transform counts for all heatmaps
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsdMatrix <- assay(vsd)

# Heatmap of the sample-to-sample distances
# see http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-sample-to-sample-distances
mat <- as.matrix(dist(t(vsdMatrix)))
mat <- heatmap_names(mat, coldata)

title <- 'Sample distance clustering (blind)'
legend_aes <- list(
  "breaks" = c(min(mat), mean(c(min(mat), max(mat))), max(mat)),
  "labels" = c("low",    "distance\n(euclidean)",     "high")
)
out_pdf <- sub(".png", ".pdf", out_png)

heatmap_plot(mat, title, heatmap_aes, legend_aes, out_pdf)

# convert pdf to png (required for MULTIQC)
pdf2png(out_pdf, out_png)


heatmap_aes$colors <- rev(heatmap_aes$colors)
# Heatmaps of various correlation matrices
for (method in list("spearman", "pearson")){
  name <- paste(method, 'correlation clustering')

  # Kendall is also available, but is computationally much more complex
  mat <- as.matrix(cor(vsdMatrix, method=method))
  mat <- heatmap_names(mat, coldata)

  # capitalize method name
  title <- paste0(toupper(substring(name, 1, 1)), substring(name, 2))
  legend_aes <- list(
    "breaks" = NA,
    "labels" = NA
  )
  out_pdf_cor <- sub("sample_distance_clustering", gsub(" ", "_", name), out_pdf)
  heatmap_plot(mat, title, heatmap_aes, legend_aes, out_pdf_cor)

  # convert pdf to png (required for MULTIQC)
  out_png_cor <- sub(".pdf", ".png", out_pdf_cor)
  pdf2png(out_pdf_cor, out_png_cor)
}
