suppressMessages({
  library(DESeq2)
  library(IHW)
})

# snakemake variables
threads        <- snakemake@threads[[1]]
log_file       <- snakemake@log[[1]]
counts_file    <- snakemake@input[[1]]
samples_file   <- snakemake@params[[1]]
contrast       <- snakemake@wildcards$contrast
params         <- snakemake@config$DE_params
mtp            <- snakemake@config$DE_params$multiple_testing_procedure
fdr            <- snakemake@config$DE_params$alpha_value
se             <- snakemake@config$DE_params$shrinkage_estimator
assembly       <- snakemake@wildcards$assembly
output_table   <- snakemake@output$table
output_ma_plot <- snakemake@output$ma_plot

# #test variables
# threads        <- 4
# log_file       <- '/mnt/passer/bank/experiments/2019-09/dirk_nextseq/results/log/deseq2/GRCh38.p12-diseaseTest_control_disease.diffexp.log'
# counts_file    <- '/mnt/passer/bank/experiments/2019-09/dirk_nextseq/results/gene_counts/GRCh38.p12-counts.tsv'
# samples_file   <- '/home/siebrenf/git/snakemake-workflows/workflows/rna_seq/samples.tsv'
# contrast       <- 'stageTest_1_2'
# mtp            <- 'ihw'
# fdr            <- 0.05
# se             <- 'apeglm'
# assembly       <- 'GRCh38.p12'
# output_table   <- '/mnt/passer/bank/experiments/2019-09/dirk_nextseq/results/deseq2/GRCh38.p12-diseaseTest_control_disease.diffexp.tsv'
# output_ma_plot <- '/mnt/passer/bank/experiments/2019-09/dirk_nextseq/results/deseq2/GRCh38.p12-diseaseTest_control_disease.ma_plot.svg'

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")


## parse the design contrast
# a contrast is always in the form 'batch+condition_group1_group2', where batch(+) is optional
# extracting batch, condition and groups
contr <- gsub('~' ,'' , contrast)
batch <- NA
if (grepl('\\+', contr)) {
  batch <- strsplit(contr, '\\+')[[1]][1]
  contr <- strsplit(contr, '\\+')[[1]][2]
}
contr <- strsplit(contr, '_')[[1]]
condition <- contr[1]
groups <- contr[-1]


## obtain coldata, the metadata input for DESeq2
samples <- read.delim(samples_file, row.names=1, na.strings = "")
# rename batch and condition (needed as DESeq's design cannot accept variables)
samples[,"condition"] <- samples[condition]
samples[,"batch"] <- ifelse(!is.na(batch), samples[batch], NA)
# filter for assembly and remove NAs
samples <- samples[samples$assembly == assembly & !is.na(samples[condition]), c("condition", "batch"), drop = F]

# keep only the samples belonging to our groups
coldata <- samples[(samples["condition"] == groups[1] | samples["condition"] == groups[2]), , drop = F]
# refactor in case we dropped (several) factor levels
coldata$condition <- factor(coldata$condition)
coldata$batch <- factor(coldata$batch)


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

cat('Constructing DESeq object. If an error is found below this line, it is between you and DESeq2 to to solve it!\n')
dds <- DESeqDataSetFromMatrix(countData = reduced_counts,
                              colData = coldata,
                              design = if (!is.na(batch)){~ batch + condition} else {~ condition})
cat('Finished constructing DESeq object.\n\n')
dds <- DESeq(dds)
cat('\n')

## Extract differentially expressed genes
DE_contrast <- resultsNames(dds)[2]

# correct for multiple testing
if(mtp=='IHW'){
  res <- results(dds, name=DE_contrast, alpha=fdr, filterFun=ihw)
} else {
  res <- results(dds, name=DE_contrast, alpha=fdr)
}

# log transform counts
resLFC <- lfcShrink(dds, coef=DE_contrast, res = res, type=se)


## Save the results
# create a table with all genes, and analysis results if available
expressed_genes <- as.data.frame(resLFC[order(resLFC$padj),])

missing_genes <- rownames(counts)[!(rownames(counts) %in% rownames(expressed_genes))]

unexpressed_genes <- as.data.frame(matrix(data = NA, ncol = ncol(expressed_genes), nrow = length(missing_genes)))
rownames(unexpressed_genes) <- missing_genes
colnames(unexpressed_genes) <- colnames(expressed_genes)
unexpressed_genes[,'baseMean'] <- 0

all_genes <- rbind(expressed_genes, unexpressed_genes)
write.table(all_genes, file=output_table, quote = F, sep = '\t', col.names=NA)

# plot of log fold change vs mean gene counts
plot_res <- resLFC[resLFC$padj <= fdr & !is.na(resLFC$padj), ]
plot_DEGs <- length(plot_res[,1])

svg(output_ma_plot)
plotMA(plot_res, ylim=c(-2,2),
       main = paste0(DE_contrast, '\n', plot_DEGs, ' DE genes (a = ', fdr, ')'))
invisible(dev.off())
