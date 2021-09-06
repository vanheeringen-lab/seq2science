suppressMessages({
    library(Matrix)
    library(Seurat)
})

kb_dir <- snakemake@input[[1]]
rds <- snakemale@output[[1]] 

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
files <- list.files(kb_dir, full.names = T)
spliced <- lapply(files, read_count_output, name = "spliced") 
names(spliced) <- list.files(kb_dir)
#Get a list of all Seurat objects
all_seu <-
  lapply(seq(spliced), function(x) {
    CreateSeuratObject(counts = spliced[[x]],
                       assay = "sf",
                       project = names(spliced[x]))
    
  })
# Create merged Seurat object
seu_merged <-
  merge(
    x = all_seu[[1]],
    y = all_seu[2:length(all_seu)],
    add.cell.ids = names(spliced),
    project = "merged"
  )

#Save workspace to RData object
save.image(file = rds)



