suppressMessages({
    library(Matrix)
    library(Seurat)
})

kb_dir <- dirname(snakemake@input$counts)
rds <- snakemake@output[[1]] 
assay <- snakemake@params$assay
sample <- snakemake@params$sample
isvelo <- snakemake@params$isvelo

#Set name for non velocity analysis
name <- "cells_x_genes"
if (assay == "ADT") {
  name <- "cells_x_features"
}

#Read counts into seurat object
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/counts_unfiltered/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <-
    readLines(file(paste0(dir, "/counts_unfiltered/", name, ge)))
  barcodes <-
    readLines(file(paste0(
      dir, "/counts_unfiltered/", name, ".barcodes.txt"
    )))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

#List files
counts <- read_count_output(kb_dir, name=name)
CreateSeuratObject(counts = counts,
                   assay = assay,
                   project = sample)
    
#Save workspace to RData object
save.image(file = rds)



