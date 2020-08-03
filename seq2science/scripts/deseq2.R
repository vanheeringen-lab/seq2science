suppressMessages({
  library(DESeq2)
  library(IHW)
  library(ggplot2)
})

# snakemake variables
threads         <- snakemake@threads[[1]]
log_file        <- snakemake@log[[1]]
counts_file     <- snakemake@input[[1]]
samples_file    <- snakemake@params$samples
replicates      <- snakemake@params$replicates
contrast        <- snakemake@wildcards$contrast
mtp             <- snakemake@config$DE_params$multiple_testing_procedure
fdr             <- snakemake@config$DE_params$alpha_value
se              <- snakemake@config$DE_params$shrinkage_estimator
assembly        <- snakemake@wildcards$assembly
salmon          <- snakemake@config$quantifier == "salmon"
counts_dir      <- snakemake@config$counts_dir
output          <- snakemake@output[[1]]

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
cat('contrast     <- "', contrast,     '"\n', sep = "")
cat('mtp          <- "', mtp,          '"\n', sep = "")
cat('fdr          <-',   fdr, '\n')
cat('se           <- "', se,           '"\n', sep = "")
cat('assembly     <- "', assembly,     '"\n', sep = "")
cat('salmon       <- "', salmon,       '"\n', sep = "")
cat('counts_dir   <- "', counts_dir,   '"\n', sep = "")
cat('output       <- "', output,       '"\n', sep = "")
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')


## parse the design contrast
# a contrast is always in the form 'batch+condition_group1_group2', where batch(+) is optional
batch <- NA
contr <- contrast
if (grepl('\\+', contrast)) {
  batch <- strsplit(contrast, '\\+')[[1]][1]
  contr <- strsplit(contrast, '\\+')[[1]][2]
}
contr <- strsplit(contr, '_')[[1]]
condition <- contr[1]
groups <- contr[-1]
rm(contr)


## obtain coldata, the metadata input for DESeq2
samples <- read.delim(samples_file, sep = "\t", na.strings = "", comment.char = "#", stringsAsFactors = F)
if ("replicate" %in% colnames(samples) & isTRUE(replicates)) {
  samples$replicate[is.na(samples$replicate)] <- as.character(samples$sample[is.na(samples$replicate)])
  samples <- subset(samples, !duplicated(replicate))
  row.names(samples) <- samples$replicate
} else {
  row.names(samples) <- samples$sample
}
coldata <- samples

# rename batch and condition (required as DESeq's design cannot accept variables)
coldata[,"condition"] <- coldata[condition]
coldata[,"batch"]     <- ifelse(!is.na(batch), coldata[batch], NA)

# filter for assembly and condition & order data for DESeq
coldata <- coldata[coldata$assembly == assembly & coldata$condition %in% c(groups[1], groups[2]), c("condition", "batch")]
coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = groups[2])
coldata$batch     <- factor(coldata$batch)


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

cat('Constructing DESeq object... \nTip: errors directly below this line are most likely DESeq2 related.\n\n')
dds <- DESeqDataSetFromMatrix(countData = reduced_counts,
                              colData = coldata,
                              design = if (!is.na(batch)){~ batch + condition} else {~ condition})
cat('\nFinished constructing DESeq object.\n\n')
dds <- DESeq(dds, parallel=parallel)
cat('\n')


## Extract differentially expressed genes
cat('batches & contrasts:\n', resultsNames(dds), '\n')
DE_contrast_names <- resultsNames(dds)
for (DE_contrast in DE_contrast_names){
  if (startsWith(DE_contrast, 'condition')){
    break
  }
}
cat('selected contrast:', DE_contrast, '\n\n')

# correct for multiple testing
if(mtp=='IHW'){
  res <- results(dds, name=DE_contrast, alpha=fdr, filterFun=ihw)
} else {
  res <- results(dds, name=DE_contrast, alpha=fdr)
}

# log transform counts
resLFC <- lfcShrink(dds, coef = DE_contrast, res = res, type=se)
cat('\n')


## Save the results
# create a table with all genes
expressed_genes <- as.data.frame(resLFC[order(resLFC$padj),])

missing_genes <- rownames(counts)[!(rownames(counts) %in% rownames(expressed_genes))]
unexpressed_genes <- as.data.frame(matrix(data = NA, ncol = ncol(expressed_genes), nrow = length(missing_genes)))
rownames(unexpressed_genes) <- missing_genes
colnames(unexpressed_genes) <- colnames(expressed_genes)
unexpressed_genes[,'baseMean'] <- 0

all_genes <- rbind(expressed_genes, unexpressed_genes)
write.table(all_genes, file=output, quote = F, sep = '\t', col.names=NA)
cat('DE genes table saved\n\n')


## determine DE genes
plot_res <- resLFC[resLFC$padj <= fdr & !is.na(resLFC$padj), ]
plot_DEGs <- length(plot_res[,1])
if(plot_DEGs == 0){
  cat("No differentially expressed genes found!\n")
  quit(save = "no" , status = 0)
} else {
  cat(plot_DEGs, "differentially expressed genes found!\n")
}

## generate additional files if DE genes are found

# generate MA plot (log fold change vs mean gene counts)
output_ma_plot <- sub(".diffexp.tsv", ".ma_plot.svg", output)
svg(output_ma_plot)
plotMA(plot_res, ylim=c(-2,2),
       main = paste0(contrast, '\n', plot_DEGs, ' DE genes (a = ', fdr, ')'))
invisible(dev.off())
cat('-MA plot saved\n')

if (is.na(batch)){
  # generate a PCA plot (for sample outlier detection)

  blind_vst <- varianceStabilizingTransformation(dds, blind = TRUE)

  g = plotPCA(blind_vst, intgroup="condition")

  output_pca_plot <- sub(".diffexp.tsv", ".pca_plot.svg", output)
  svg(output_pca_plot)
  plot(g + ggtitle("blind PCA") + aes(color=condition) + theme(legend.position="bottom"))
  invisible(dev.off())
  cat('-PCA plot saved\n')

} else {
  # generate a PCA plots before and after batch correction

  blind_vst <- varianceStabilizingTransformation(dds, blind = TRUE)
  nonblind_vst <- varianceStabilizingTransformation(dds, blind = FALSE)
  # model the effect of batch correction in the correlation heatap. DESeq2 applies this internally as well.
  mat <- assay(nonblind_vst)
  mat <- limma::removeBatchEffect(mat, nonblind_vst$batch)
  batchcorr_vst <- nonblind_vst
  assay(batchcorr_vst) <- mat

  g1 = plotPCA(blind_vst, intgroup=c("condition", "batch"))
  g2 = plotPCA(batchcorr_vst, intgroup=c("condition", "batch"))

  output_pca_plots <- sub(".diffexp.tsv", ".pca_plot_%01d.svg", output)
  svg(output_pca_plots)
  if(length(levels(blind_vst$batch)) < 7){
    plot(g1 + ggtitle("blind PCA - color by condition") + aes(color=condition, shape=batch) + theme(legend.position="bottom"))
    plot(g1 + ggtitle("blind PCA - color by batch") + aes(color=batch) + theme(legend.position="bottom"))

    plot(g2 + ggtitle("batch corrected PCA - color by condition") + aes(color=condition, shape=batch) + theme(legend.position="bottom"))
    plot(g2 + ggtitle("batch corrected PCA - color by batch") + aes(color=batch) + theme(legend.position="bottom"))
  } else {
    plot(g1 + ggtitle("blind PCA - color by condition") + aes(color=condition) + theme(legend.position="bottom"))
    plot(g1 + ggtitle("blind PCA - color by batch") + aes(color=batch) + theme(legend.position="bottom"))

    plot(g2 + ggtitle("batch corrected PCA - color by condition") + aes(color=condition) + theme(legend.position="bottom"))
    plot(g2 + ggtitle("batch corrected PCA - color by batch") + aes(color=batch) + theme(legend.position="bottom"))
  }
  invisible(dev.off())
  cat('-PCA plots saved\n\n')

  # Generate the batch corrected counts
  # (for downstream tools that do not model batch effects)
  batch_corrected_counts <- function(dds) {
    #' accepts a large DESeqDataSet, normalizes and removes the batch effects, and returns a batch corrected count matrix

    nonblind_vst <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
    mat <- SummarizedExperiment::assay(nonblind_vst)
    mat <- limma::removeBatchEffect(mat, nonblind_vst$batch)

    # mat contains the normalized and batch corrected variant stabilized data (mat = vst.fn(batch_corr_counts))
    # next, we invert this function (from DESeq2::getVarianceStabilizedData) to get normalized and batch corrected counts (batch_corr_counts)

    # vst.fn <- function(q) {
    #     log((1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] *
    #         q + 2 * sqrt(coefs["asymptDisp"] * q * (1 + coefs["extraPois"] +
    #         coefs["asymptDisp"] * q)))/(4 * coefs["asymptDisp"]))/log(2)
    # }
    coefs <- attr(DESeq2::dispersionFunction(dds), "coefficients")
    x <- mat * log(2)
    x <- exp(1)^x * (4 * coefs["asymptDisp"])
    x <- (x - (1 + coefs["extraPois"]))/2
    batch_corr_counts <- x^2 / (coefs["asymptDisp"] * (coefs["extraPois"] + 2 * x + 1))
    return(batch_corr_counts)
  }

  batch_corr_counts <- batch_corrected_counts(dds)
  output_batch_corr_counts <- sub(".diffexp.tsv", ".batch_corr_counts.tsv", output)
  write.table(batch_corr_counts, file=output_batch_corr_counts, quote = F, sep = '\t', col.names=NA)
  cat('-batch corrected counts saved\n')

  # if quantified with salmon, generate the batch corrected TPMs as well
  if(salmon){
    RDS_path <- file.path(counts_dir, paste0(assembly, "-se.rds"))
    se <- readRDS(RDS_path)
    sg <- se$sg
    gene_lengths <- assay(sg, "length")

    counts_to_tpm <- function(mat.counts, mat.gene_lengths) {
      #' snippet from https://support.bioconductor.org/p/91218/ to convert counts to TPM

      x <- mat.counts / mat.gene_lengths
      mat.tpm <- t( t(x) * 1e6 / colSums(x) )
      return(mat.tpm)
    }
    batch_corr_tpm <- counts_to_tpm(batch_corr_counts, gene_lengths)

    output_batch_corr_tpm <- sub(".diffexp.tsv", ".batch_corr_tpm.tsv", output)
    write.table(batch_corr_tpm, file=output_batch_corr_tpm, quote = F, sep = '\t', col.names=NA)
    cat('-batch corrected TPMs saved\n')
  }
}
