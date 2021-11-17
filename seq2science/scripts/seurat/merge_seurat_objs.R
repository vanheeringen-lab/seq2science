library(Seurat)

#Load parameters
seu_dir <- paste0(snakemake@config$result_dir,"/seurat/kallistobus")
rds <- snakemake@output[[1]] 
#isvelo <- snakemake@params$isvelo
#iskite <- snakemake@params$iskite

dir <- normalizePath(seu_dir, mustWork = TRUE)
output_folders <- list.files(dir, 
                             recursive = TRUE, include.dirs = FALSE, full.names = TRUE )

#Read all RDS object from directory
res <- sapply(output_folders, readRDS)

#Get the names for each object
cell.ids <- c()
for (i in 1:length(res)){
  id <- res[[i]]@project.name
  cell.ids <- c(cell.ids, id)
}
names(res) <- cell.ids

#Merge Seurat object and keep cell id.
if (length(res) > 1) {
seu_obj <-
  merge(x = res[[1]],
        y = res[2:length(res)],
        add.cell.ids = names(res),
        project = "merge")
}  else {
  seu_obj <- res[[1]]
}

saveRDS(seu_obj, file = rds)
