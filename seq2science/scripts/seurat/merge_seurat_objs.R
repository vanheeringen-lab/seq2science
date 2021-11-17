library(Seurat)

#Snakemake variables
seu_dir <- paste0(snakemake@config$result_dir,"/seurat/kallistobus")
rds <- snakemake@output[[1]] 
log_file <- snakemake@log[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('seu_dir      <- "', seu_dir,         '"\n', sep = "")
cat('rds          <- "', rds,             '"\n', sep = "")
cat('\n')

#isvelo <- snakemake@params$isvelo
#iskite <- snakemake@params$iskite

# List RDS files
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

#Merge Seurat object and use cell id as prefix.
if (length(res) > 1) {
seu_obj <-
  merge(x = res[[1]],
        y = res[2:length(res)],
        add.cell.ids = names(res),
        project = "merge")
}  else {
  seu_obj <- res[[1]]
}

#Save merged object in RDS format
saveRDS(seu_obj, file = rds)
