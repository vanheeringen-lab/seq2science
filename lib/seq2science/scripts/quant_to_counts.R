# source: https://github.com/csoneson/ARMOR/blob/master/scripts/run_tximeta.R

suppressMessages({
  library(tidyverse)
  library(SingleCellExperiment)
})

# snakemake variables
linked_txome <- snakemake@input$linked_txome
samples      <- snakemake@input$cts
sample_names <- snakemake@params$names
assembly     <- snakemake@wildcards$assembly
out_counts   <- snakemake@output$counts
out_tpms     <- snakemake@output$tpms
out_lengths  <- snakemake@output$lengths
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
cat('sample_names <- "', sample_names, '"\n', sep = "")
cat('assembly     <- "', assembly,     '"\n', sep = "")
cat('out_counts   <- "', out_counts,   '"\n', sep = "")
cat('out_tpms     <- "', out_tpms,     '"\n', sep = "")
cat('out_lengths  <- "', out_lengths,  '"\n', sep = "")
cat('out_SCE      <- "', out_SCE,      '"\n', sep = "")
cat('log_file     <- "', log_file,     '"\n', sep = "")
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')


## Load linked_txome.json
tximeta::loadLinkedTxome(linked_txome)

coldata <- data.frame(files = file.path(samples, 'quant.sf'), names = sample_names, stringsAsFactors = F, check.names = F)

## import annotated abundances in transcript level
st <- tximeta::tximeta(
    coldata=coldata,
    txOut=TRUE,                # output transcripts (default)
    # skipMeta=TRUE,           # meta = required for transcript outputs
    # skipSeqinfo=TRUE,        # lookup sizes
    useHub=FALSE,              # lookup similar indexes
    # markDuplicateTxps=TRUE,  # mark and track
    cleanDuplicateTxps=TRUE,   # fix
)

## Summarize to gene level
sg <- tximeta::summarizeToGene(st)


## This section deviates from the source script
## It outputs non-normalized matrixes

## Save TPM matrix
TPM <- data.frame(assay(sg, "abundance"), stringsAsFactors = F, check.names = F) %>% rownames_to_column("gene")
write.table(TPM, file=out_tpms, quote = F, sep = '\t', row.names = F)

## Save gene counts matrix
counts <- assay(sg, "counts")
mode(counts) <- "integer"
counts <- data.frame(counts, stringsAsFactors = F, check.names = F) %>% rownames_to_column("gene")
write.table(counts, file=out_counts, quote = F, sep = '\t', row.names = F)

## Save gene length matrix
lengths <- data.frame(assay(sg, "length"), stringsAsFactors = F, check.names = F) %>% rownames_to_column("gene")
write.table(lengths, file=out_lengths, quote = F, sep = '\t', row.names = F)

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
  data.frame()

## Change the row names in sg to have geneID__geneSymbol
rownames(sg) <- paste(rowData(sg)$gene_id, rowData(sg)$gene_name, sep = "__")

## Coerce the object from SummarizedExperiment to SingleCellExperiment
st <- as(st, "SingleCellExperiment")
sg <- as(sg, "SingleCellExperiment")

## Save single cell experiment object
saveRDS(list(st = st, sg = sg), file = out_SCE)
