"assess the fraction of artifactual reads to normal read duplication (due to natural over-sequencing of highly expressed genes)."

suppressMessages({
  library(dupRadar)
})

# snakemake variables
log_file       <- snakemake@log[[1]]
bam_file       <- snakemake@input$bam
gtf_file       <- snakemake@input$gtf
strandedness   <- snakemake@params$strandedness
paired         <- snakemake@params$paired
sample         <- snakemake@wildcards$sample
threads        <- snakemake@threads[[1]]
out_plot       <- snakemake@output[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('log_file     <- "', log_file,      '"\n', sep = "")
cat('bam_file     <- "', bam_file,      '"\n', sep = "")
cat('gtf_file     <- "', gtf_file,      '"\n', sep = "")
cat('strandedness <- "', strandedness,  '"\n', sep = "")
cat('paired       <- "', paired,        '"\n', sep = "")
cat('sample       <- "', sample,        '"\n', sep = "")
cat('threads      <- ',  threads,        '\n', sep = "")
cat('out_plot     <- "', out_plot,      '"\n', sep = "")
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')

# # '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
# if (stranded == "reverse"){
#   stranded <- 2
# } else if (stranded == "forward"){
#   stranded <- 1
# } else {
#   stranded <- 0
# }

# .libPaths("/home/siebrenf/miniconda3/envs/s2s/envs/dupradar/lib/R/library")
# The call parameters:
# bam <- "./mm10-external-Caro-colon-RNA-13-18710.samtools-coordinate.bam"
# gtf <- "./mm10.annotation.gtf"
# stranded <- 2       # '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
# paired   <- TRUE
# threads  <- 4
# outdir <- "."
# name <- gsub(".samtools-coordinate.bam","",basename(bam))
# attach(dupRadar_examples)  # example dm
# out_plot <- file.path(outdir, paste0(name, "_dupRadar.png"))

# Analysis
dm <- analyzeDuprates(bam_file, gtf_file, strandedness, paired, threads, verbose = TRUE)

# Plot
png(file=out_png, width=3000,  height=1000)
par(mfrow=c(1,3), oma=c(0, 0, 2, 0))
duprateExpDensPlot(dm)
title("Density plot")
expressionHist(dm)
title("Boxplot")
duprateExpBoxplot(dm)
title("Histplot")
mtext(sample, outer = TRUE, cex = 1.5)
dev.off()
