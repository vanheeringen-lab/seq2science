suppressMessages({
    library(Matrix)
    library(Seurat)
})

kb_dir <- dirname(snakemake@input$counts)
rds <- snakemake@output[[1]] 
assay <- snakemake@params$assay
sample <- snakemake@params$sample
isvelo <- snakemake@params$isvelo
iskite <- snakemake@params$iskite

#Set name for non velocity analysis

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
# Create Seurat objects based on input workflow
if (iskite) {
  assign(paste("seu.",sample,".ADT"),CreateSeuratObject(counts = read_count_output(kb_dir, name="cells_x_features"), assay = "ADT", project = sample))
} else if (isvelo) {
   assign(paste("seu.",sample,".sf"), CreateSeuratObject(counts = read_count_output(kb_dir, name="spliced"), assay = "sf", project = sample))                   
   assign(paste("seu.",sample,".uf"), CreateSeuratObject(counts = read_count_output(kb_dir, name="unspliced"), assay = "uf", project = sample))
} else {
  assign(paste("seu.",sample,".RNA"), CreateSeuratObject(counts = read_count_output(kb_dir, name="cells_x_genes"), assay = "RNA", project = sample))
}

#Save workspace to RData object
save.image(file = rds)



