# source: https://github.com/csoneson/ARMOR/blob/master/scripts/run_tximeta.R

suppressMessages({
  library(tidyverse)
  library(tximeta)
  library(readr)
  library(SingleCellExperiment)
})

# snakemake variables
linked_txome <- snakemake@input$linked_txome
samples      <- snakemake@input$cts
assembly     <- snakemake@wildcards$assembly
out_matrix   <- snakemake@output$counts
out_SCE      <- snakemake@output$SCE
log_file     <- snakemake@log[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('linked_txome <- "', linked_txome, '"\n', sep = "")
cat('samples      <- "', samples,      '"\n', sep = "")
cat('assembly     <- "', assembly,     '"\n', sep = "")
cat('out_matrix   <- "', out_matrix,   '"\n', sep = "")
cat('out_SCE      <- "', out_SCE,      '"\n', sep = "")
cat('log_file     <- "', log_file,     '"\n', sep = "")
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')


## Load linked_txome.json
loadLinkedTxome(linked_txome)

samplenames <- gsub(paste0(assembly, '-'), '', basename(samples))
coldata <- data.frame(names = samplenames, files = file.path(samples, 'quant.sf'), stringsAsFactors = F, check.names = F)

## import annotated abundances in transcript level
st <- tximeta::tximeta(coldata)

## Summarize to gene level
sg <- summarizeToGene(st)


## This section deviates from the source script
## It outputs a un-normalized counts matrix

## Extract a count table
counts <- assay(sg, "counts")

## Convert to int
mode(counts) <- "integer"

## Save gene counts matrix
counts <- data.frame(counts, stringsAsFactors = F, check.names = F) %>% rownames_to_column("gene")
write.table(counts, file=out_matrix, quote = F, sep = '\t', row.names = F)

## Returning to source script


## If rowData(st)$gene_id is a CharacterList, convert it to character to allow 
## the joining below
if (is(rowData(st)$gene_id, "CharacterList")) {
  if (any(vapply(rowData(st)$gene_id, length, 1) > 1)) {
    warning("Some elements of rowData(st)$gene_id consisted of more than one",
            "object. Only the first one is retained.")
  }
  rowData(st)$gene_id <- vapply(rowData(st)$gene_id, function(w) w[[1]], "")
}

## If rowData(st)$tx_id is of class integer, replace it with the tx_name 
## column
if (is(rowData(st)$tx_id, "integer")) {
  rowData(st)$tx_id <- rowData(st)$tx_name
}

## Add gene information, e.g. gene_name, entrezid, ... (if provided) to
## transcript-level SE
rowData(st) <- rowData(st) %>%
  data.frame() %>%
  dplyr::left_join(data.frame(rowData(sg))) %>%
  DataFrame()

## Change the row names in sg to have geneID__geneSymbol
rownames(sg) <- paste(rowData(sg)$gene_id, rowData(sg)$gene_name, sep = "__")

## Coerce the object from SummarizedExperiment to SingleCellExperiment
st <- as(st, "SingleCellExperiment")
sg <- as(sg, "SingleCellExperiment")

## Save single cell experiment object
saveRDS(list(st = st, sg = sg), file = out_SCE)
