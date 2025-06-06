#!/usr/bin/env Rscript

# input method
if (exists("snakemake")){
  # parse seq2science input
  scripts_dir        <- snakemake@params$scripts_dir
  deseq_init         <- file.path(scripts_dir, "run_as_rule.R")
  output             <- snakemake@output$diffexp
  output_ma_plot     <- snakemake@output$maplot
  output_vol_plot    <- snakemake@output$volcanoplot
  output_pca_plot    <- snakemake@output$pcaplot
} else {
  # parse command line input
  cli_args           <- commandArgs(trailingOnly=F)
  this_script        <- sub("--file=", "", cli_args[grep("--file=", cli_args)])
  scripts_dir        <- dirname(this_script)
  deseq_init         <- file.path(scripts_dir, "run_as_standalone.R")
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
coldata               <- coldata[c("condition", "batch")]

# log involved samples
cat('samples per group:\n')
cat('  - target    (', condition, groups[1],'): ', rownames(coldata[coldata$condition %in% groups[1],]),'\n')
cat('  - reference (', condition, groups[2],'): ', rownames(coldata[coldata$condition %in% groups[2],]),'\n\n')

# determine if we need to run batch correction on the whole assembly
output_batch_corr_counts <- sub(paste0(contrast, ".diffexp.tsv"), paste0(batch, ".batch_corr_counts.tsv"), output, fixed=TRUE)
output_batch_corr_pca    <- sub(paste0(contrast, ".diffexp.tsv"), paste0(batch, ".batch_corr_pca.png"), output, fixed=TRUE)
output_batch_corr_tpm    <- sub(paste0(contrast, ".diffexp.tsv"), paste0(batch, ".batch_corr_tpm.tsv"), output, fixed=TRUE)

# filter unused samples & order data for DESeq
coldata <- coldata[!is.na(coldata$condition),]
cat('Using all (', length(rownames(coldata[!is.na(coldata$condition), ])) ,') ', sep = "")
cat('samples labelled in the condition column ("', condition, '") ', sep = "")
cat('to calculate the dispersions:', rownames(coldata[!is.na(coldata$condition), ]), '\n\n')
if (!is.na(batch)) {
  if (any(is.na(coldata$batch))) {
    cat('Error: all samples labelled in the condition column ("', condition, '") ', sep = "")
    cat('need a label in the batch column ("', batch,'").\n', sep = "")
    cat('Unlabelled samples:', rownames(coldata[is.na(coldata$batch), ]), '\n\n')
    quit(save = "no" , status = 1)
  }
}
coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = groups[2])
coldata$batch     <- factor(coldata$batch)


## filter counts to speed up DESeq
counts <- read.table(counts_file, row.names = 1, header = T, stringsAsFactors = F, sep = '\t', quote= "", check.names = F)
reduced_counts <- counts[rowSums(counts) > 0, rownames(coldata)]

## DESeq2
design <- if (!is.na(batch)){~ batch + condition} else {~ condition}
dds <- run_deseq2(reduced_counts, coldata, design, threads, single_cell)


## Extract differentially expressed genes
cat('batches & contrasts:\n', DESeq2::resultsNames(dds), '\n\n')
DE_contrast <- paste("condition", groups[1], "vs", groups[2], sep="_")
cat('selected contrast:', DE_contrast, '\n\n')

# correct for multiple testing
if (mtp=='IHW') {
  res <- DESeq2::results(dds, name=DE_contrast, alpha=fdr, filterFun=ihw)
} else {
  res <- DESeq2::results(dds, name=DE_contrast, alpha=fdr)
}

# log transform counts
resLFC <- DESeq2::lfcShrink(dds, coef=DE_contrast, res=res, type=se)
cat('\n')


## Save the results
# save a diffexp table with all genes, not just the expressed genes
save_complete_diffexp(resLFC, counts, output)
cat('DE genes table saved\n\n')


## determine DE genes
n_DEGs <- length(resLFC[resLFC$padj <= fdr & !is.na(resLFC$padj), ][,1])
if (n_DEGs == 0) {
  cat("No differentially expressed genes found!\n")
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

png(output_ma_plot, width = 250, height = 250, units='mm', res = 300)
DESeq2::plotMA(
  resLFC,
  alpha = fdr,
  ylab = 'log2 fold change',
  main = title
)
invisible(dev.off())
cat('-MA plot saved\n')


png(output_vol_plot, width = 250, height = 250, units='mm', res = 300)
EnhancedVolcano::EnhancedVolcano(
  resLFC,
  lab = rownames(resLFC),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = title
)
invisible(dev.off())
cat('-volcano plot saved\n')

if (is.na(batch)) {
  # generate a PCA plot (for sample outlier detection)

  blind_vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  g <- DESeq2::plotPCA(blind_vst, intgroup="condition")

  png(output_pca_plot, width = 465, height = 225, units='mm', res = 300)
  plot(g + ggtitle("blind PCA") + aes(color=condition) + theme(legend.position="bottom"))
  invisible(dev.off())
  cat('-PCA plot saved\n')

} else {
  # generate a PCA plots before and after batch correction (all samples and contrast samples)

  blind_vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  batchcorr_vst <- batch_corrected_vst(dds)

  # batch corrected PCA (all samples)
  if (!file.exists(output_batch_corr_pca)) {
    g1 <- DESeq2::plotPCA(blind_vst, intgroup=c("condition", "batch"))
    g2 <- DESeq2::plotPCA(batchcorr_vst, intgroup=c("condition", "batch"))

    # prevent plot() from generating Rplot.pdf in the working directory
    pdf(NULL)

    # color by batch/condition. up to 6 shapes can be displayed too.
    condition_aes <- if (length(levels(blind_vst$batch)) < 7) {aes(color=condition, shape=batch)} else {aes(color=condition)}
    batch_aes <- if (length(levels(blind_vst$condition)) < 7) {aes(color=batch, shape=condition)} else {aes(color=batch)}

    plot1 <- plot(g1 + ggtitle("blind PCA - color by condition") + condition_aes + theme(legend.position="bottom"))
    plot2 <- plot(g1 + ggtitle("blind PCA - color by batch") + batch_aes + theme(legend.position="bottom"))

    plot3 <- plot(g2 + ggtitle("batch corrected PCA - color by condition") + condition_aes + theme(legend.position="bottom"))
    plot4 <- plot(g2 + ggtitle("batch corrected PCA - color by batch") + batch_aes + theme(legend.position="bottom"))

    png(output_batch_corr_pca, width = 465, height = 225, units='mm', res = 300)
    gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2)
    invisible(dev.off())
    cat('-batch corrected PCA plots saved\n')
  } else {
    cat('-batch corrected PCA plots already exists\n')
  }

  # batch corrected PCA (contrast samples)
  # subset batch corrected data to contrast samples
  blind_vst <- blind_vst[,rownames(coldata)[coldata$condition %in% c(groups[1], groups[2])]]
  batchcorr_vst <- batchcorr_vst[,rownames(coldata)[coldata$condition %in% c(groups[1], groups[2])]]

  g1 <- DESeq2::plotPCA(blind_vst, intgroup=c("condition", "batch"))
  g2 <- DESeq2::plotPCA(batchcorr_vst, intgroup=c("condition", "batch"))

  # prevent plot() from generating Rplot.pdf in the working directory
  pdf(NULL)

  # color by batch/condition. up to 6 shapes can be displayed too.
  condition_aes <- if (length(levels(blind_vst$batch)) < 7) {aes(color=condition, shape=batch)} else {aes(color=condition)}
  batch_aes <- aes(color=batch, shape=condition)

  plot1 <- plot(g1 + ggtitle("blind PCA - color by condition") + condition_aes + theme(legend.position="bottom"))
  plot2 <- plot(g1 + ggtitle("blind PCA - color by batch") + batch_aes + theme(legend.position="bottom"))

  plot3 <- plot(g2 + ggtitle("batch corrected PCA - color by condition") + condition_aes + theme(legend.position="bottom"))
  plot4 <- plot(g2 + ggtitle("batch corrected PCA - color by batch") + batch_aes + theme(legend.position="bottom"))

  png(output_pca_plot, width = 465, height = 225, units='mm', res = 300)
  gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2)
  invisible(dev.off())
  cat('-PCA plots saved\n\n')

  # batch corrected counts (all samples)
  # (for downstream tools that do not model batch effects)
  if (!file.exists(output_batch_corr_counts)) {
    batch_corr_counts <- batch_corrected_counts(dds)

    bcc <- cbind(gene = rownames(batch_corr_counts), batch_corr_counts)  # consistent with other tables
    write.table(bcc, file=output_batch_corr_counts, quote = F, sep = '\t', row.names = F)
    cat('-batch corrected counts saved\n')
  } else {
    cat('-batch corrected counts already exists\n')
  }

  # batch corrected TPM (all samples)
  # (if quantified with salmon)
  if (salmon & !file.exists(output_batch_corr_tpm)) {
    if (!exists("batch_corr_counts")) {
      batch_corr_counts <- read.table(output_batch_corr_counts, row.names = 1, header = T, stringsAsFactors = F, sep = '\t', check.names = F)
    }
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
