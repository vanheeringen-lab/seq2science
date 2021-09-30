#!/usr/bin/env Rscript

# input method
if (exists("snakemake")){
  # parse seq2science input
  scripts_dir     <- snakemake@params$scripts_dir
  deseq_init      <- file.path(scripts_dir, "run_as_rule.R")
  output          <- snakemake@output$diffexp
  output_ma_plot  <- snakemake@output$maplot
  output_pca_plot <- snakemake@output$pcaplot
} else {
  # parse command line input
  cli_args        <- commandArgs(trailingOnly=F)
  this_script     <- sub("--file=", "", cli_args[grep("--file=", cli_args)])
  scripts_dir     <- dirname(this_script)
  deseq_init      <- file.path(scripts_dir, "run_as_standalone.R")
  output_ma_plot  <- sub(".diffexp.tsv", ".ma_plot.pdf", output)
  output_pca_plot <- sub(".diffexp.tsv", ".pca_plot.pdf", output)
}
deseq_utils <- file.path(scripts_dir, "utils.R")
source(deseq_init)
source(deseq_utils)

# parse the design contrast
ret <- parse_contrast(contrast)
batch <- ret$batch          # a column in samples or NA
condition <- ret$condition  # a column in samples
groups <- ret$groups        # >1 field in samples[condition]


## obtain coldata, the metadata input for DESeq2
samples <- parse_samples(samples_file, assembly, replicates)

# rename batch and condition (required as DESeq's design cannot accept variables)
coldata <- samples
coldata[,"condition"] <- coldata[condition]
coldata[,"batch"]     <- if (!is.na(batch)) { coldata[batch] } else { NA }
coldata <- coldata[c("condition", "batch")]

# determine if we need to run batch correction on the whole assembly
output_batch_corr_counts <- sub(paste0(contrast, ".diffexp.tsv"), paste0(batch, "+", condition, ".batch_corr_counts.tsv"), output, fixed=TRUE)
output_batch_corr_tpm <- sub(paste0(contrast, ".diffexp.tsv"), paste0(batch, "+", condition, ".batch_corr_tpm.tsv"), output, fixed=TRUE)
no_batch_correction_required <- is.na(batch) | (file.exists(output_batch_corr_counts) & (!salmon | file.exists(output_batch_corr_tpm)))

# filter out unused conditions & order data for DESeq
if (no_batch_correction_required) {
  coldata <- coldata[coldata$condition %in% c(groups[1], groups[2]),]
} else {
  cat('\nbatch correction dataset selected\n\n')
  # for batch corrected counts we want all samples marked in the batch column
  coldata <- coldata[!is.na(coldata$batch),]
}
coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = groups[2])
coldata$batch     <- factor(coldata$batch)


## filter counts to speed up DESeq
counts <- read.table(counts_file, row.names = 1, header = T, stringsAsFactors = F, sep = '\t', check.names = F)
reduced_counts <- counts[rowSums(counts) > 0, colnames(counts) %in% rownames(coldata)]


## DESeq2
design <- if (!is.na(batch)){~ batch + condition} else {~ condition}
dds <- run_deseq2(reduced_counts, coldata, design, threads)


## Extract differentially expressed genes
cat('batches & contrasts:\n', resultsNames(dds), '\n\n')
DE_contrast <- paste("condition", groups[1], "vs", groups[2], sep="_")
cat('selected contrast:', DE_contrast, '\n\n')

# correct for multiple testing
if (mtp=='IHW') {
  res <- results(dds, name=DE_contrast, alpha=fdr, filterFun=ihw)
} else {
  res <- results(dds, name=DE_contrast, alpha=fdr)
}

# log transform counts
resLFC <- lfcShrink(dds, coef = DE_contrast, res = res, type=se)
cat('\n')


## Save the results
# save a diffexp table with all genes, not just the expressed genes
save_complete_diffexp(resLFC, counts, output)
cat('DE genes table saved\n\n')


## determine DE genes
n_DEGs <- length(resLFC[resLFC$padj <= fdr & !is.na(resLFC$padj), ][,1])
if (n_DEGs == 0) {
  cat("No differentially expressed genes found!\n")
  quit(save = "no" , status = 0)
} else {
  cat(n_DEGs, "differentially expressed genes found!\n")
}

## generate additional files if DE genes are found

# generate MA plot (log fold change vs mean gene counts)
b <- ifelse(is.na(batch), '', ', batch corrected')
title <- paste0(
  groups[1], ' vs ', groups[2], '\n',
  n_DEGs, ' of ', nrow(reduced_counts), ' DE (a = ', fdr, b, ')'
)
png(output_ma_plot)
DESeq2::plotMA(
  resLFC,
  alpha = fdr,
  ylab = 'log2 fold change',
  main = title
)
invisible(dev.off())
cat('-MA plot saved\n')

if (is.na(batch)) {
  # generate a PCA plot (for sample outlier detection)

  blind_vst <- varianceStabilizingTransformation(dds, blind = TRUE)
  g <- DESeq2::plotPCA(blind_vst, intgroup="condition")

  png(output_pca_plot)
  plot(g + ggtitle("blind PCA") + aes(color=condition) + theme(legend.position="bottom"))
  invisible(dev.off())
  cat('-PCA plot saved\n')

} else {
  # generate a PCA plots before and after batch correction

  blind_vst <- varianceStabilizingTransformation(dds, blind = TRUE)
  g1 <- DESeq2::plotPCA(blind_vst, intgroup=c("condition", "batch"))

  batchcorr_vst <- batch_corrected_vst(dds)
  g2 <- DESeq2::plotPCA(batchcorr_vst, intgroup=c("condition", "batch"))

  # color by batch/condition. up to 6 shapes can be displayed too.
  condition_aes <- if (length(levels(blind_vst$batch)) < 7) {aes(color=condition, shape=batch)} else {aes(color=condition)}
  batch_aes <- if (length(levels(blind_vst$condition)) < 7) {aes(color=batch, shape=condition)} else {aes(color=batch)}

  output_pca_plots <- sub(".pca_plot.png", ".pca_plot_%01d.pdf", output_pca_plot)
  png(output_pca_plots)
  plot(g1 + ggtitle("blind PCA - color by condition") + condition_aes + theme(legend.position="bottom"))
  plot(g1 + ggtitle("blind PCA - color by batch") + batch_aes + theme(legend.position="bottom"))

  plot(g2 + ggtitle("batch corrected PCA - color by condition") + condition_aes + theme(legend.position="bottom"))
  plot(g2 + ggtitle("batch corrected PCA - color by batch") + batch_aes + theme(legend.position="bottom"))
  invisible(dev.off())
  cat('-PCA plots saved\n\n')

  # Generate the batch corrected counts
  # (for downstream tools that do not model batch effects)
  if (!file.exists(output_batch_corr_counts)) {
    batch_corr_counts <- batch_corrected_counts(dds)

    bcc <- cbind(gene = rownames(batch_corr_counts), batch_corr_counts)  # consistent with other tables
    write.table(bcc, file=output_batch_corr_counts, quote = F, sep = '\t', row.names = F)
    cat('-batch corrected counts saved\n')
  } else {
    cat('-batch corrected counts already exists\n')
  }

  # if quantified with salmon, generate the batch corrected TPMs as well
  if (salmon & !file.exists(output_batch_corr_tpm)) {
    lengths_file <- sub("-counts.tsv", "-gene_lengths.tsv", counts_file)
    gene_lengths <- read.table(lengths_file, row.names = 1, header = T, stringsAsFactors = F, sep = '\t', check.names = F)
    gene_lengths <- gene_lengths[rownames(batch_corr_counts),colnames(batch_corr_counts)]

    batch_corr_tpm <- counts2tpm(batch_corr_counts, gene_lengths)
    write.table(batch_corr_tpm, file=output_batch_corr_tpm, quote = F, sep = '\t', col.names=NA)
    cat('-batch corrected TPMs saved\n')
  } else if (salmon & file.exists(output_batch_corr_tpm)) {
    cat('-batch corrected TPMs already exists\n')
  }
}
