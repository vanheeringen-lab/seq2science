suppressMessages({
  library(dupRadar)
})

# snakemake variables
log_file       <- snakemake@log[[1]]
bam_file       <- snakemake@input$bam
gtf_file       <- snakemake@input$gtf
strandedness   <- as.numeric(snakemake@params$strandedness)
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
cat('strandedness <- ',  strandedness,   '\n', sep = "")
cat('paired       <- ',  paired,         '\n', sep = "")
cat('sample       <- "', sample,        '"\n', sep = "")
cat('threads      <- ',  threads,        '\n', sep = "")
cat('out_plot     <- "', out_plot,      '"\n', sep = "")
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')

# Add examples
good_example = file.path(dirname(out_plot), "good_example.png")
bad_example = file.path(dirname(out_plot), "bad_example.png")
if (!file.exists(good_example)){
  attach(dupRadar_examples)

  png(file=good_example, width=1344,  height=672)
  par(mfrow=c(1,2), oma=c(0, 0, 2, 0))
  duprateExpDensPlot(dm)
  title("Density plot")
  # expressionHist(dm)
  # title("Boxplot")
  duprateExpBoxplot(dm)
  title("Histplot")
  mtext("good example", outer = TRUE, cex = 1.5)
  dev.off()

  png(file=bad_example, width=1344,  height=672)
  par(mfrow=c(1,2), oma=c(0, 0, 2, 0))
  duprateExpDensPlot(dm.bad)
  title("Density plot")
  # expressionHist(dm.bad)
  # title("Boxplot")
  duprateExpBoxplot(dm.bad)
  title("Histplot")
  mtext("bad example", outer = TRUE, cex = 1.5)
  dev.off()

  rm(dm)
}

# Analysis
dm <- analyzeDuprates(bam_file, gtf_file, strandedness, paired, threads, verbose = TRUE)

# Plot
png(file=out_plot, width=1344,  height=672)  # same dimensions as their vignette
par(mfrow=c(1,2), oma=c(0, 0, 2, 0))
duprateExpDensPlot(dm)
title("Density plot")
# expressionHist(dm)
# title("Boxplot")
duprateExpBoxplot(dm)
title("Histplot")
mtext(sample, outer = TRUE, cex = 1.5)
dev.off()
