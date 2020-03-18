suppressMessages({
  library(R.utils)
  library(tximeta)
})

# snakemake variables
salmonindex <- snakemake@input$index_dir
source      <- snakemake@params$source
organism    <- snakemake@params$organism
release     <- snakemake@params$release
genome      <- snakemake@wildcards$assembly
fasta       <- snakemake@input$fasta
gtf         <- snakemake@input$gtf
output      <- snakemake@output$index
log_file    <- snakemake@log[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('salmonindex <- "', salmonindex, '"\n', sep = "")
cat('source      <- "', source,      '"\n', sep = "")
cat('organism    <- "', organism,    '"\n', sep = "")
cat('release     <- ',  release,     '\n',  sep = "")
cat('genome      <- "', genome,      '"\n', sep = "")
cat('fasta       <- "', fasta,       '"\n', sep = "")
cat('gtf         <- "', gtf,         '"\n', sep = "")
cat('output      <- "', output,      '"\n', sep = "")
cat('log_file    <- "', log_file,    '"\n', sep = "")
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')


## create a symlink to the gtf with Ensembl naming scheme (required for tximeta)
fake_gtf_path <- file.path(dirname(output), paste0(organism, '.', genome, '.', release, '.gtf'))
cat('Renaming GTF:\n')
createLink(fake_gtf_path, gtf)
cat('\n')

## Creating linked transcriptome
cat('Creating linked transcriptome:\n')
makeLinkedTxome(indexDir = salmonindex,
                source = source,
                organism = organism,
                release = release,
                genome = genome,
                fasta = fasta,
                gtf = fake_gtf_path, # instead of the real gtf
                write = TRUE,
                jsonFile = output)
