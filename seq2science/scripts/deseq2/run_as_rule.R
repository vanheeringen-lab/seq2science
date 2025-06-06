suppressMessages({
  library(ggalt)  # for the volcano plots
})

# snakemake variables
threads         <- snakemake@threads[[1]]
log_file        <- snakemake@log[[1]]
counts_file     <- snakemake@input[[1]]
samples_file    <- snakemake@params$samples
replicates      <- snakemake@params$replicates
contrast        <- snakemake@wildcards$contrast
assembly        <- snakemake@wildcards$assembly
mtp             <- snakemake@config$deseq2$multiple_testing_procedure
fdr             <- snakemake@config$deseq2$alpha_value
se              <- snakemake@config$deseq2$shrinkage_estimator
single_cell     <- snakemake@config$deseq2$single_cell
salmon          <- snakemake@params$salmon
output          <- snakemake@output[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('threads      <- ',  threads,      '\n',  sep = "")
cat('log_file     <- "', log_file,     '"\n', sep = "")
cat('counts_file  <- "', counts_file,  '"\n', sep = "")
cat('samples_file <- "', samples_file, '"\n', sep = "")
cat('replicates   <- ',  replicates,   '\n',  sep = "")
cat('contrast     <- "', contrast,     '"\n', sep = "")
cat('assembly     <- "', assembly,     '"\n', sep = "")
cat('mtp          <- "', mtp,          '"\n', sep = "")
cat('fdr          <- ',  fdr,          '\n',  sep = "")
cat('se           <- "', se,           '"\n', sep = "")
cat('single_cell  <- ',  single_cell,  '\n',  sep = "")
cat('salmon       <- ',  salmon,       '\n',  sep = "")
cat('output       <- "', output,       '"\n', sep = "")
cat('\n')

cat('sessioninfo:\n')
sessionInfo()
cat('\n')
