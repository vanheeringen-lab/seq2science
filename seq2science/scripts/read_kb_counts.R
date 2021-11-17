suppressMessages({
    library(Matrix)
    library(Seurat)
})

#Set name for non velocity analysis
prep_cell_meta <- function(seu, sample_sheet) {
  blacklist <- c("result_folder","descriptive_name","technical_replicates")
  meta <- colnames(sample_sheet[setdiff(names(sample_sheet),blacklist)])
  sample.meta <- sample_sheet[sample_sheet$sample==seu@project.name,]
  sample.meta <- sample.meta[meta]
  sample.meta <- sample.meta[rep(seq_len(nrow(sample.meta)), each = ncol(seu)),]
  rownames(sample.meta) <- rownames(FetchData(seu,"ident"))
  sample.meta <- sample.meta[match(rownames(FetchData(seu,"ident")),rownames(sample.meta)),]
  return(sample.meta)
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
#Load parameters
kb_dir <- dirname(snakemake@input$counts)
rds <- snakemake@output[[1]] 
assay <- snakemake@params$assay
sample <- snakemake@params$sample
isvelo <- snakemake@params$isvelo
iskite <- snakemake@params$iskite

sample_sheet <- tryCatch(
    read.table(snakemake@config$samples, sep = '\t', header = TRUE)
)

# Create Seurat objects based on input workflow
if (iskite) {
  seu <- CreateSeuratObject(counts = read_count_output(kb_dir, name="cells_x_features"), assay = "ADT", project = sample)
  seu@meta.data <- prep_cell_meta(seu, sample_sheet)
  saveRDS(seu, file = rds)

} else if (isvelo) {
   seu.sf <- CreateSeuratObject(counts = read_count_output(kb_dir, name="spliced"), assay = "sf", project = sample)
   seu.uf <- CreateSeuratObject(counts = read_count_output(kb_dir, name="unspliced"), assay = "uf", project = sample)
   seu.sf@meta.data <- prep_cell_meta(seu.sf, sample_sheet)
   seu.uf@meta.data <- prep_cell_meta(seu.uf, sample_sheet)
   saveRDS(c(seu.sf, seu.uf), file = rds) 
    
} else {
  seu <- CreateSeuratObject(counts = read_count_output(kb_dir, name="cells_x_genes"), assay = "RNA", project = sample)
  seu@meta.data <- prep_cell_meta(seu, sample_sheet)
  saveRDS(seu, file = rds)

}





