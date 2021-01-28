read_kb_counts <- function(dir, name, barcode_file, remove_bc=TRUE) {
  # Loading scRNA-seq count matrix
  #
  # Generates a combined scRNA-seq matrix from the ouput of 
  # the Kallisto | Bustools wrapper. Needs the location of
  # the count tables to combine. Returned matrix consists of 
  # genes (with gene names) in the rows, and cells (with 
  # barcode identifier) in the columns. 
  #
  # dir: location of the "_output" folders generated 
  #   by the kb-wrapper
  # name: the name of count assay you want to load (e.g. 
  #   spliced, unspliced or cells_x_genes, the latter when 
  #   running without velocity)

  # barcode_file is a 2-column file with a well-id its respective DNA barcode
  # t2g_file is the transcript-to-gene file created with the kb-wrapper genome index

  # Change only the parameters in this block:
  #----------------------------#
  # Location of needed files
  barcode_file <- barcode_file
  #----------------------------#

  ## Loading packages & files ##
  library(Matrix)
  library(tidyr)
  library(dplyr)
  
  ## Generate matrix ##
  
  # Iterate over all _output folders, generating matrix per plate
  # combining matrices by gene matches
  dir <- normalizePath(dir, mustWork = TRUE)
  output_folders <- list.files(dir, 
                               recursive = FALSE, include.dirs = TRUE)
  i <- 1
  for (folder in output_folders){
    print(paste("Reading:",folder))
    plate <- paste0(dir, "/", folder, "/counts_unfiltered/", name)
    m <- readMM(paste0(plate, ".mtx"))
    m <- Matrix::t(m)
    m <- as(m, "dgCMatrix")
    # the matrix has genes in columns, cells in rows, 
    # stored in .genes.txt and .barcodes.txt
    genes <- as.vector(read.table(file(paste0(plate, ".genes.txt")))[,1])
    barcodes <- as.vector(read.table(file(paste0(plate, ".barcodes.txt")))[,1])
    # retrieve unique plate-id from folder
    #platename <- gsub("GRCh38-", "", folder)
    colnames(m) <- paste(barcodes, folder, sep = "_")
    rownames(m) <- genes
    # create a combined matrix for all plates in the folder
    if (i == 1) {
      combined <- m
    } else if (identical(rownames(combined),rownames(m)) == TRUE){
      # Only binds the matrices if genes are identical and in the same order
      combined <- cbind(combined, m)
    }
    i <- i + 1
  }
  if (remove_bc){
    ## Replace cell barcodes for well identifier ##
    # barcode file contains the well identifier and corresponding DNA barcode
    plate_order <- read.table(barcode_file, sep = "\t", col.names = c("well","barcode"))
    # generate a data.frame to match barcode and wellid 
    cells <- data.frame("cell" = colnames(combined))
    cells$barcode <- gsub("_.*", "", cells$cell)
    cells$well <- plate_order$well[match(cells$barcode, plate_order$barcode)]
    # Remove DNA barcode and add wellid 
    cells$cell_id <- paste(gsub("^.*?_", "", cells$cell), cells$well, sep = "_")
    cells$cell_id <- gsub("-", "_", cells$cell_id)
    # replace cell names of the count matrix
    colnames(combined) <- cells$cell_id
  }
  return(combined)
}