args <- commandArgs(trailingOnly=TRUE)

# help message
if ("--help" %in% args || "-h" %in% args) {
  cat('\nUSAGE:\n  ./deseq2.R [1]contrast [2]/path/to/samples_file [3]/path/to/counts_file [4]/path/to/out_dir  \n\n')
  quit(save = "no" , status = 0)
}

# conda environment
required_packages <- c("DESeq2", "BiocParallel", "IHW", "ggplot2")
installed_packages <- installed.packages()[,"Package"]
for (pkg in required_packages){
  if ( !(pkg %in% installed_packages) ){
    cat('Please activate a conda environment with DESeq2, BiocParallel, IHW and ggplot2.  \n\n')
    quit(save = "no" , status = 0)
  }
}

# parse arguments
if (!length(args) %in% c(4,5,6)) {
  cat("Four arguments expected: contrast, samples_file, counts_file and outdir. Optional: assembly, single_cell.")
  quit(save = "no" , status = 0)
}
contrast     <- args[1]
samples_file <- args[2]
counts_file  <- args[3]
out_dir      <- args[4]
assembly     <- args[5]
single_cell  <- ifelse(length(args)==6, as.logical(args[6]), FALSE)

if (!length(strsplit(contrast, '_')[[1]]) >= 3){
  cat("Contrast must contain a column name and two fields separated with an underscore (_).  \n\n")
  quit(save = "no" , status = 0)
}

if (!file.exists(samples_file)){
  cat("Cannot find the samples_file.  \n\n")
  quit(save = "no" , status = 0)
}

if (!file.exists(counts_file)){
  cat("Cannot find the counts_file.  \n\n")
  quit(save = "no" , status = 0)
}

if (!dir.exists(out_dir)){
  dir.create(out_dir, showWarnings = F, recursive = T)
}

# variables required in the core script
samples            <- read.delim(samples_file, sep = "\t", na.strings = "", comment.char = "#", stringsAsFactors = F, row.names = "sample", check.names = F)
replicates         <- "technical_replicates" %in% colnames(samples)                          # always merge replicates if "technical_replicates" exists
assembly           <- ifelse(assembly %in% samples$assembly, assembly, samples$assembly[1])  # use the first assembly
mtp                <- "BH"                                                                   # \
fdr                <- 0.1                                                                    #  |-default options only
se                 <- "apeglm"                                                               # /
salmon             <- FALSE                                                                  # only work with counts data
threads            <- 4
output             <- file.path(out_dir, paste0(assembly, "-", contrast, ".diffexp.tsv"))
output_ma_plot     <- sub(".diffexp.tsv", ".ma_plot.png", output)
output_vol_plot    <- sub(".diffexp.tsv", ".volcano_plot.png", output)
output_pca_plot    <- sub(".diffexp.tsv", ".pca_plot.png", output)

# load libraries
suppressMessages({
  library(ggalt)  # for the volcano plots
})
