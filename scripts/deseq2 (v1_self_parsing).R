# suppressMessages({
  library(DESeq2)
  library(IHW)
#   library(pheatmap)
#   library(RColorBrewer)
#   library(tidyverse)
# })

# # snakemake variables
# threads        <- snakemake@threads[[1]]
# log_file       <- snakemake@log[[1]]
# counts_file    <- snakemake@input[[1]]
# samples_file   <- snakemake@params[[1]]
# contrast       <- snakemake@wildcards$contrast
# mtp            <- snakemake@config$$diffexp$deseq2$mtp
# fdr            <- snakemake@config$diffexp$deseq2$alpha
# se             <- snakemake@config$diffexp$deseq2$se
# assembly       <- snakemake@wildcards$assembly
# output_table   <- snakemake@output$table
# output_ma_plot <- snakemake@output$ma_plot

#test variables
threads        <- 4
log_file       <- '/mnt/passer/bank/experiments/2019-09/dirk_nextseq/results/log/deseq2/GRCh38.p12-diseaseTest_control_disease.diffexp.log'
counts_file    <- '/mnt/passer/bank/experiments/2019-09/dirk_nextseq/results/gene_counts/GRCh38.p12-counts.tsv'
samples_file   <- '/home/siebrenf/git/snakemake-workflows/workflows/rna_seq/samples.tsv'
contrast       <- '~ stageTest_1_all'
mtp            <- 'ihw'
fdr            <- 0.05
se             <- 'apeglm'
assembly       <- 'GRCh38.p12'
output_table   <- '/mnt/passer/bank/experiments/2019-09/dirk_nextseq/results/deseq2/GRCh38.p12-diseaseTest_control_disease.diffexp.tsv'
output_ma_plot <- '/mnt/passer/bank/experiments/2019-09/dirk_nextseq/results/deseq2/GRCh38.p12-diseaseTest_control_disease.ma_plot.svg'

# # log all console output
# log <- file(log_file, open="wt")
# sink(log)
# sink(log, type="message")



## parse the design contrast
# A contrast is always in the form '~ batch + condition_group1_group2', where batch and groups are optional.
contr <- gsub(' |~' ,'' , contrast)
if (grepl('\\+', contr)) {
  batch <- strsplit(contr, '\\+')[[1]][1]
  contr <- strsplit(contr, '\\+')[[1]][2]
} else {
  batch <- NA
}
contr <- strsplit(contr, '_')[[1]]
condition <- contr[1]
groups <- contr[-1] #ifelse(length(contr) == 1, NA, contr[2:length(contr)])

# safety check, contrast parsing
test <- try(all(
  !is.null(condition),
  length(groups) <= 2
), silent = T)
if (!(test == T)) stop('design contrast is messed up yo!')



## compare & combine contrast with samples.tsv
# samples.tsv has columns named after the batch and condition. The condition column has rows containing the groups.
samples = read.delim(samples_file, row.names=1, na.strings = "")

# safety check, presence of batch and conditions in samples
if (!(condition %in% colnames(samples))) stop(paste0('Error! Condition ', condition, ' not found in samples'))
if (!is.na(batch)) {if (!(batch %in% colnames(samples))) stop(paste0('Error! Batch ', batch, ' not found in samples'))}


## obtain coldata, the metadata input for DESeq2
# rename batch and presence (needed as DESeq's design cannot accept variables)
samples[,"condition"] <- samples[condition]
samples$condition <- factor(samples$condition)
samples[,"batch"] <- ifelse(!is.na(batch), samples[batch], NA)
samples$batch <- factor(samples$batch)

samples = samples[samples$assembly == assembly & !is.na(samples[condition]), c("condition", "batch"), drop = F]
if (is.na(groups[1])) {
  #A vs B, only two groups exist
  coldata = samples[,, drop = F]
  reflvl = levels(samples$condition)[1]
  
  # safety check, expecting only 2 groups
  if(length(levels(coldata$condition)) > 2) stop(paste0('Your design contrast only contains a column name (', condition,'), but the column contains more than two groups (', length(levels(coldata$condition)), ' found)!'))

} else if (length(groups) == 2 & !('all' %in% groups)) {
  #A vs B, ignore C, D, etc.
  coldata = samples[(samples["condition"] == groups[1] | samples["condition"] == groups[2]), , drop = F]
  reflvl = groups[1]
  
} else if (length(groups) == 1 | 'all' %in% groups) {
  #A vs all.
  coldata = samples[,, drop = F]
  reflvl = ifelse(groups[1] != 'all', groups[1], groups[2])

} else {
  #safety check. this should never happen...
  stop(paste0('Error in design condition. Expecting output in the form Header_group1_group2, found ', paste(contrast, sep = '_', collapse='_'), '. Read the wiki for more information on DE designs.'))
}
coldata$condition = relevel(samples$condition, reflvl)



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

dds <- DESeqDataSetFromMatrix(countData = reduced_counts,
                              colData = coldata,
                              design = if (!is.na(batch)){~ batch + condition} else {~ condition})
dds <- DESeq(dds)


## Extract differentially expressed genes
# most experiments will have one comparisson: group1 vs group2
# for these cases resultsNames(dds) will be length 2, thus the loop will run only once.
# however, if the design contrast is group1 vs all, there will be (length(all)-1) loops.
for (DE_contrast in resultsNames(dds)[2:length(resultsNames(dds))]){
  # correct for multiple testing
  if(mtp=='ihw'){
    res <- results(dds, name=DE_contrast, alpha=fdr, filterFun=ihw)
  } else {
    res <- results(dds, name=DE_contrast, alpha=fdr)
  }

  # log transform counts
  resLFC <- lfcShrink(dds, coef=DE_contrast, res = res, type=se)

  # create a table with all genes and analysis results if available
  expressed_genes <- as.data.frame(resLFC[order(resLFC$padj),])
  
  newrows = rownames(counts)[!(rownames(counts) %in% rownames(expressed_genes))]
  unexpressed_genes <- as.data.frame(matrix(data = NA, ncol = ncol(expressed_genes), nrow = length(newrows)))
  rownames(unexpressed_genes) <- newrows
  colnames(unexpressed_genes) <- colnames(expressed_genes)
  unexpressed_genes[,'baseMean'] <- 0
  
  all_genes <- rbind(expressed_genes, unexpressed_genes)
  
  # plot of log fold change vs mean gene counts
  plot_res <- resLFC[resLFC$padj <= fdr & !is.na(resLFC$padj), ]
  plot_DEGs = length(plot_res[,1])
  
  # the table and plot are ready and need to be saved. If the design contrast contains multiple experiments, these need to be split up.
  if (length(resultsNames(dds)) == 2){
    # save the table
    write.table(all_genes, file=output_table, quote = F, sep = '\t', col.names=NA)
    
    # save the plot
    svg(output_ma_plot)
    plotMA(plot_res, ylim=c(-2,2),
           main = paste0(DE_contrast, '\n', plot_DEGs, ' DE genes (a = ', fdr, ')'))
    dev.off()
  } else {
    if (!file.exists(output_table)){
      # create dummy files (to satisfy snakemake)
      #dir.create(dname)
      dirmessage='This contrast design contains multiple conditions which can be found in the subdirectory.'
      write.table(dirmessage, file=output_table, quote = F, sep = '\t', col.names=NA)
      
      svg(output_ma_plot)
      plot.new()
      dev.off()
    }
    dname = dirname(output_table)
    
    # save the table
    fname = paste0(dname, '/', assembly, '-', DE_contrast, '.diffexp.tsv')
    write.table(all_genes, file=fname, quote = F, sep = '\t', col.names=NA)
    
    # save the plot
    pname = paste0(dname, '/', assembly, '-', DE_contrast, '.ma_plot.svg')
    svg(pname)
    plotMA(plot_res, ylim=c(-2,2),
           main = paste0(DE_contrast, '\n', DEGs, ' DE genes (a = ', fdr, ')'))
    dev.off()
  }
}



# # correct for multiple testing
# contrast = resultsNames(dds)[2]
# if (mtp=='ihw'){
#   res <- results(dds, name=contrast, filterFun=ihw, alpha=fdr)
# } else {
#   res <- results(dds, name=contrast, alpha=fdr)
# }
# # log transform counts
# resLFC <- lfcShrink(dds, coef=contrast, res = res, type="apeglm")
# 
# # save all genes with DE analysis data if possible
# expressed_genes <- as.data.frame(resLFC[order(resLFC$padj),])
# 
# newrows = rownames(counts)[!(rownames(counts) %in% rownames(resLFC))]
# unexpressed_genes <- as.data.frame(matrix(data = NA, ncol = ncol(expressed_genes), nrow = length(newrows)))
# rownames(unexpressed_genes) <- newrows
# colnames(unexpressed_genes) <- colnames(expressed_genes)
# unexpressed_genes[,'baseMean'] <- 0
# 
# all_genes <- rbind(expressed_genes, unexpressed_genes)
# 
# write.table(all_genes, file=output_table, quote = F, sep = '\t', col.names=NA)
# 
# # plot significant DE genes
# plot_res <- resLFC[resLFC$padj <= fdr & !is.na(resLFC$padj), ]
# DEGs = length(plot_res[,1])
# svg(output_ma_plot)
# plotMA(plot_res, ylim=c(-2,2),
#        main = paste0(contrast, '\n', DEGs, ' DE genes (a = ', fdr, ')'))
# dev.off()










# # library(dendextend)
# # library(pheatmap)
# # library(DESeq2)
# # library("BiocParallel")
# # library("IHW")
# # library(tidyverse)
# #
# # input
# # file='/bank/experiments/2019-09/dirk_nextseq/results/gene_counts/GRCh38.p12-counts.tsv'
# # meta='/bank/experiments/2019-09/dirk_nextseq/Samples_Dirk.csv'
# # 
# # # output
# # bias_clust_path='/bank/experiments/2019-09/dirk_nextseq/results/deseq/GRCh38.p12-bias_detection_clustering.pdf'
# # clust_path='/bank/experiments/2019-09/dirk_nextseq/results/deseq/GRCh38.p12-sample_clustering.pdf'
# # save_path='/bank/experiments/2019-09/dirk_nextseq/results/deseq/GRCh38.p12-'
# 
# 
# 
# # Load the counts table, sort by sample and remove empty rows
# counts = read.table(input_file, row.names = 1, header = T, stringsAsFactors = F, sep = '\t', check.names = F)
# 
# counts = t(counts)
# counts = counts[order(row.names(counts)), colSums(counts) > 0]
# counts = t(counts)
# 
# # Load the metadata, sort by sample name, keep only samples present in counts
# metadata = read.table(meta, row.names = 1, header = T, stringsAsFactors = F) 
# 
# metadata = metadata[order(colnames(counts)), ]
# 
# 
# 
# # ### hierarchical clustering
# # 
# # ### raw data - can only be used to determine biases or big mistakes (switched samples etc.)
# # # convert the counts to a dendrogram
# # library(dendextend)
# # dend = t(counts) %>%
# #   dist() %>% # several methods available
# #   hclust() %>% # several methods available
# #   as.dendrogram()
# # 
# # # add factor levels to metadata (used for coloring)
# # DF <- as.data.frame(unclass(metadata), row.names = rownames(metadata))
# # 
# # pdf(file=bias_clust_path)
# # # color by flowcell, run, category, expecting_CPF and granules
# # for (col in colnames(metadata)){
# #   # determine color schemes
# #   l = c()
# #   for (label in labels(dend)){
# #     l2 = as.numeric(DF[label,col])
# #     l =c(l, l2)
# #   }
# #   labels_colors(dend) <- l
# # 
# #   # plot
# #   par(mar = c(2,1,1,5))
# #   plot(dend, horiz = TRUE, main = paste0('color by ', col))
# #   legend("topleft", legend = unique(as.character(DF[,col])), fill = unique(as.numeric(DF[,col])), cex = 0.7)
# # }
# # dev.off()
# 
# ### corrected count values - used to inspect differences or similarities between samples
# library(DESeq2)
# library("BiocParallel")
# register(MulticoreParam(10))
# library("RColorBrewer")
# library(pheatmap)
# 
# coldata = metadata[, 'category', drop=FALSE]
# colnames(coldata) <- c("condition")
# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = coldata,
#                               design = ~ 1) #dummy design, we only want the sample correlations
# dds <- DESeq(dds)
# 
# # correct for library size imbalance, and transform the counts & plot
# # https://www.biostars.org/p/351551/
# # http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
# if (nrow(coldata) > 30){
#   main = 'Sample distances (vst)'
#   log_counts=varianceStabilizingTransformation(dds)
# } else {
#   main = 'Sample distances (rlog)'
#   log_counts <- rlog(dds, blind = TRUE)
# }
# log_counts=assay(log_counts)
# 
# ## print
# pdf(file=clust_path, width = 11)
# ## cluster based on correlation
# # sampleCorMatrix <- cor(log_counts)
# # colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)
# # pheatmap(sampleCorMatrix,
# #          cex=1, main = 'Sample correlation',
# #          annotation_col=coldata,
# #          col=colors)
# 
# ## cluster based on poisson distance
# # library("PoiClaClu")
# # poisd <- PoissonDistance(t(counts(dds)))
# # samplePoisDistMatrix <- as.matrix( poisd$dd )
# # rownames(samplePoisDistMatrix) <- colnames(counts(dds))
# # colnames(samplePoisDistMatrix) <- colnames(counts(dds))
# # pheatmap(samplePoisDistMatrix,
# #          cex=1, main = main,
# #          annotation_col=coldata,
# #          clustering_distance_rows = poisd$dd,
# #          clustering_distance_cols = poisd$dd,
# #          col = colors)
# 
# ## cluster based on sample distance
# sampleDistMatrix <- as.matrix(dist(t(log_counts)))
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          cex=1, main = main,
#          annotation_col=coldata,
#          col=colors)
# dev.off()
# 
# 
# 
# ### DE analysis
# library("IHW")
# 
# # add the conditions requested
# metadata$cond1 = ifelse(grepl("control", metadata$category), "control", ifelse(grepl("active.*diseas", metadata$category), 'active_disease', NA))
# metadata$cond2 = ifelse(grepl("remission", metadata$category), "remission", ifelse(grepl("active.*diseas", metadata$category), 'active_disease', NA))
# metadata$cond3 = ifelse(grepl("remission", metadata$category), "remission", ifelse(grepl("(disease).*control", metadata$category), "disease_control", NA))
# 
# metadata$cond4 = ifelse("no" == metadata$expecting_CPF, "no", ifelse("yes" == metadata$expecting_CPF, "yes", NA))
# metadata$cond5 = ifelse("no_or_less" == metadata$expecting_CPF, "no_or_less", ifelse("yes" == metadata$expecting_CPF, "yes", NA))
# metadata$cond6 = ifelse("no" == metadata$expecting_CPF, "no", ifelse("no_or_less" == metadata$expecting_CPF, "no_or_less", NA))
# 
# metadata$cond7 = ifelse("no" == metadata$granules, "no", 'yes')
# 
# metadata = metadata[,grepl('cond.', colnames(metadata))]
# 
# 
# 
# # loop over all conditions
# for (condition in colnames(metadata)){
#   # coldata tells DESeq which condition to analyse
#   coldata = metadata[,condition, drop=FALSE]
#   colnames(coldata) <- c("condition")
#   coldata = coldata[rowSums(is.na(coldata)) != ncol(coldata), , drop=FALSE]
#   coldata$condition = factor(coldata$condition)
#   
#   # only take samples involved in the condition
#   reduced_counts = counts[, colnames(counts) %in% rownames(coldata)]
#   reduced_counts = reduced_counts[rowSums(reduced_counts) > 0,]
#   
#   # DESeq black box
#   dds <- DESeqDataSetFromMatrix(countData = reduced_counts,
#                                 colData = coldata,
#                                 design = ~ condition)
#   dds <- DESeq(dds)
#   
#   ### Indentify DE genes 
#   # both the multiple testing and shrinkage method need to be cited (see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#plot-counts)
#   
#   # accepted false discovery rate:
#   alpha_value = 0.05
#   
#   contrast = resultsNames(dds)[2]
#   # multiple-testing procedure: Independent Hypothesis Weighting (instead of Benjamini-Hochberg, as this takes the condition into account)
#   res <- results(dds, name=contrast, filterFun=ihw, alpha=alpha_value)
#   # log transform using apeglm (may be switched for ashr)
#   resLFC <- lfcShrink(dds, coef=contrast, res = res, type="apeglm")
#   DEGs = sum(resLFC$padj < alpha_value, na.rm=TRUE)
#   
#   # plot DE genes
#   plot_path = paste0(save_path, "DE_", condition, '.pdf')
#   pdf(file=plot_path)
#   plotMA(resLFC, ylim=c(-2,2), alpha = alpha_value,
#          main = paste0(contrast, '\n', DEGs, ' DE genes (a = 0.05)'))
#   dev.off()
#   
#   # save all DE genes to file with their log2 fold-changes and adjusted p-values
#   table_path = paste0(save_path, "DE_", condition, '.tsv')
#   DEG_list = as.matrix(resLFC[resLFC$padj < alpha_value & !is.na(resLFC$padj), ])
#   DEG_list = DEG_list[order(DEG_list[,2]),]
#   write.table(DEG_list, table_path, quote = F, sep = '\t', col.names=NA)
# }
# dev.off()
