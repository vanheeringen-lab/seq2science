suppressMessages({
    library(Matrix)
    library(Seurat)
    library(dplyr)
})

kb_dir <- snakemake@input$counts
rds <- snakemake@output[[1]] 
assay <- snakemake@params$assay
sample <- snakemake@params$sample
name <- "cell_x_genes"

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



