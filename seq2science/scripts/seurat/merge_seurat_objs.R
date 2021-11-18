suppressMessages({
    library(Seurat)
})

# Snakemake variables
seu_dir <- paste0(snakemake@config$result_dir,"/seurat/kallistobus")
rds <- snakemake@output[[1]] 
isvelo <- snakemake@params$isvelo
log_file <- snakemake@log[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('seu_dir      <- "', seu_dir,         '"\n', sep = "")
cat('rds          <- "', rds,             '"\n', sep = "")
cat('isvelo       <- "', isvelo,             '"\n', sep = "")

cat('\n')

# List RData files
dir <- normalizePath(seu_dir, mustWork = TRUE)
output_folders <- list.files(dir, 
                             recursive = TRUE, include.dirs = FALSE, full.names = TRUE, pattern = "^.*_seu_obj\\.RData$" )
                             
# Helper function to merge Seurat objects
merge_seu_objects <- function(res){
  cell.ids <- lapply(res,attr, which="project.name")
  names(res) <- unlist(cell.ids)
  if (length(res) > 1) {
    seu_obj <-
      merge(x = res[[1]],
            y = res[2:length(res)],
            add.cell.ids = names(res),
            project = "merge")
  }  else {
    seu_obj <- res[[1]]
  }
}                             

# Read velocity data if detected (Lamanno workflow argument)
if (isvelo) {
  spliced <- c()
  unspliced <- c()
  for (i in 1:length(output_folders)){
    data <- readRDS(output_folders[i])
    spliced <- c(spliced,data["sf"])
    unspliced <- c(unspliced, data["uf"])  
  }
  # Merge spliced counts
  seu_obj_sf <- merge_seu_objects(spliced)
  # Merged unspliced counts
  seu_obj_uf <- merge_seu_objects(unspliced)
  # Save workspace
  seu_velo_merged <- c(seu_obj_sf, seu_obj_uf)
  names(seu_velo_merged) <- c("sf_merged","uf_merged")
  saveRDS(seu_velo_merged, file = rds)
} else {
  # Read all RDS object from directory
  res <- sapply(output_folders, readRDS)
  # Merge counts
  seu_obj <- merge_seu_objects(res)
  # Save merged object in RDS format
  saveRDS(seu_obj, file = rds)
}
