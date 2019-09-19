suppressMessages({
  library(tximeta)
  library(R.utils)
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

# create a symlink to the gtf with Ensembl naming scheme (required for tximeta)
fake_gtf_path <- file.path(dirname(output), paste0(organism, '.', genome, '.', release, '.gtf'))
createLink(fake_gtf_path, gtf, overwrite = F)

# create the linked transcriptome
makeLinkedTxome(indexDir = salmonindex,
                source = source,
                organism = organism,
                release = release,
                genome = genome,
                fasta = fasta,
                gtf = fake_gtf_path, # instead of the real gtf
                write = TRUE, 
                jsonFile = output)
