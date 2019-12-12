suppressMessages({
  library(R.utils)
  library(tximeta)
})

## snakemake variables
salmonindex <- snakemake@input$index_dir
source      <- snakemake@params$source
organism    <- snakemake@params$organism
release     <- snakemake@params$release
genome      <- snakemake@wildcards$assembly
fasta       <- snakemake@input$fasta
gtf         <- snakemake@input$gtf
output      <- snakemake@output$index
log_file    <- snakemake@log[[1]]

## log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")


cat('Sessioninfo:\n')
sessionInfo()
cat('\n')

## create a symlink to the gtf with Ensembl naming scheme (required for tximeta)
fake_gtf_path <- file.path(dirname(output), paste0(organism, '.', genome, '.', release, '.gtf'))
cat('Renaming GTF:\n')
createLink(fake_gtf_path, gtf)
cat('\n')

#
#
#
#temp: test input untill tximeta is updated on conda
# https://support.bioconductor.org/p/126276/
# https://github.com/mikelove/tximeta/issues/21
salmonindex <- "/bank/genomes/GRCh38.p13/index/salmon"
source <- "Ensembl"
organism <- "tximeta_file"
release <- "42"
genome <-"GRCh38.p13"
fasta <- "/bank/genomes/GRCh38.p13/GRCh38.p13.transcripts.fa"
fake_gtf_path <- "/bank/genomes/GRCh38.p13/index/tximeta/tximeta_file.GRCh38.p13.42.gtf"
output <- "/bank/genomes/GRCh38.p13/index/tximeta/linked_txome.json"

gtf = '/bank/genomes/GRCh38.p13/GRCh38.p13.gtf'
createLink(fake_gtf_path, gtf)
# end of temp
#
#
#

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
