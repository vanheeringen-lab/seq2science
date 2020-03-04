suppressMessages({
  library(DESeq2)
  library(IHW)
})

# snakemake variables
threads         <- snakemake@threads[[1]]
log_file        <- snakemake@log[[1]]
counts_file     <- snakemake@input[[1]]
samples_file    <- snakemake@params[[1]]
contrast        <- snakemake@wildcards$contrast
mtp             <- snakemake@config$DE_params$multiple_testing_procedure
fdr             <- snakemake@config$DE_params$alpha_value
se              <- snakemake@config$DE_params$shrinkage_estimator
assembly        <- snakemake@wildcards$assembly
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
cat('contrast     <- "', contrast,     '"\n', sep = "")
cat('mtp          <- "', mtp,          '"\n', sep = "")
cat('fdr          <-',   fdr, '\n')
cat('se           <- "', se,           '"\n', sep = "")
cat('assembly     <- "', assembly,     '"\n', sep = "")
cat('output       <- "', output,       '"\n', sep = "")
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
coldata <- read.delim(samples_file, row.names=1, na.strings = "")

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


## determine DE genes and plot if found
plot_res <- resLFC[resLFC$padj <= fdr & !is.na(resLFC$padj), ]
plot_DEGs <- length(plot_res[,1])
if(plot_DEGs == 0){
  cat("No differentially expressed genes found! Skipping plot generation...\n")
  quit(save = "no" , status = 0)
} else {
  cat(plot_DEGs, "differentially expressed genes found! Plotting MA and PCA...\n")
}

# generate MA plot (log fold change vs mean gene counts)
output_ma_plot <- sub(".diffexp.tsv", ".ma_plot.svg", output)
svg(output_ma_plot)
plotMA(plot_res, ylim=c(-2,2),
       main = paste0(contrast, '\n', plot_DEGs, ' DE genes (a = ', fdr, ')'))
invisible(dev.off())
cat('-MA plot saved\n')

# transform the data and generate a PCA plot (for outlier detection)
log_counts <- vst(dds, blind = TRUE)

output_pca_plot <- sub(".diffexp.tsv", ".pca_plot.svg", output)
svg(output_pca_plot)
plotPCA(log_counts, intgroup="condition")
invisible(dev.off())
cat('-PCA plot saved\n')
