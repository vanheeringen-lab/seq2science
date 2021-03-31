suppressMessages({
  library(edgeR)
})

# snakemake variables
log_file         <- snakemake@log[[1]]
counts_tsv       <- snakemake@input[[1]]
method           <- snakemake@wildcards$normalisation
norm_counts_tsv  <- snakemake@output[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# read count table
counts <- read.delim(counts_tsv, sep = "\t", na.strings = "", comment.char = "#", stringsAsFactors = F, , row.names=1)

# normalize
tryCatch(
    expr = {
        dgelist <- calcNormFactors(DGEList(counts), method=method)
        norm_counts <- cpm(dgelist)
    },
    error = function(e){
        norm_counts <<- data.frame(counts)
        norm_counts[] <<- "NA"
        print(paste("Something went wrong when converting the count table with method", method, "."))
    }
)

# and save
cat(paste("# The number of reads under each peak, cpm", method, "normalized\n", sep=" "), file=norm_counts_tsv)
write.table(norm_counts, file=norm_counts_tsv, quote=FALSE, sep='\t', col.names=NA, append=TRUE)
