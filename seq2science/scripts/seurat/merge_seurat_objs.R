suppressMessages({
    library(Seurat)
})

#Snakemake variables
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

# List RDS files
dir <- normalizePath(seu_dir, mustWork = TRUE)
output_folders <- list.files(dir, 
                             recursive = TRUE, include.dirs = FALSE, full.names = TRUE )
                             
#Helper function to merge seurat objects
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

#Read velocity data if detected
if(isvelo) {
  res <- c()
  unspliced <- c()
  for (i in 1:length(output_folders)){
    data <- readRDS(output_folders[i])
    res <- c(res,data[[1]])
    unspliced <- c(unspliced, data[[2]])  
  }
  #Merge spliced counts
  seu_obj <- merge_seu_objects(res)
  #Merged unspliced counts
  seu_obj_uns <- merge_seu_objects(unspliced)
  #Save workspace
  saveRDS(c(seu_obj,seu_obj_uns), file = rds)
} else {
  #Read all RDS object from directory
  res <- sapply(output_folders, readRDS)
  # Merge counts
  seu_obj <- merge_seu_objects(res)
  #Save merged object in RDS format
  saveRDS(seu_obj, file = rds)
}
