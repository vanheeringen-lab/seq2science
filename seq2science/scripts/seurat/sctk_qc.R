# This script is a lightweight adaption of the SCTK_runQC.R script for seq2science
# https://github.com/compbiomed/singleCellTK/blob/master/exec/SCTK_runQC.R

suppressMessages({
  library(Seurat)
  library(singleCellTK)
  library("BiocParallel")
})

# Snakemake variables
rds_in <- snakemake@input$rds_raw
out_dir <-snakemake@params$outdir
log_file <- snakemake@log[[1]]
sample <-  snakemake@params$sample
isvelo <- snakemake@params$isvelo
replicates <- snakemake@params$replicates
data_type <- snakemake@config$sctk$data_type
mito_set <- snakemake@config$sctk$mito_set
detect_cell <- snakemake@config$sctk$detect_cell
detect_mito <- snakemake@config$sctk$detect_mito
cell_calling <- snakemake@config$sctk$cell_calling
rds_out <- file.path(out_dir, "seu_obj_sctk.RData",    fsep="/" )
qc_out <-  file.path(out_dir, "SCTK_CellQC_summary.csv",    fsep="/" )
pdf_out <- file.path(out_dir, "SCTK_DropletQC_figures.pdf", fsep="/" )
numCores <- 4

# Log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# Log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('log_file         <- "', log_file,         '"\n', sep = "")
cat('sample           <- "', sample,           '"\n', sep = "")
cat('replicates       <- "', replicates,       '"\n', sep = "")
cat('isvelo           <- "', isvelo,           '"\n', sep = "")
cat('rds_in           <- "', rds_in,           '"\n', sep = "")
cat('out_dir          <- "', out_dir,          '"\n', sep = "")
cat('rds_out          <- "', out_dir,          '"\n', sep = "")
cat('qc_out           <- "', qc_out,           '"\n', sep = "")
cat('pdf_out          <- "', pdf_out,          '"\n', sep = "")
cat('data_type        <- "', data_type,        '"\n', sep = "")
cat('detect_mito      <- "', detect_mito,      '"\n', sep = "")
cat('mito_set         <- "', mito_set,         '"\n', sep = "")
cat('detect_cell      <- "', detect_cell,      '"\n', sep = "")
cat('cell_calling     <- "', cell_calling,     '"\n', sep = "")


# Setup parallel type
# https://github.com/compbiomed/singleCellTK/blob/master/exec/SCTK_runQC.R
parallelType <- "MulticoreParam"
Params <- list()

if (numCores > 1) {
    if (numCores > parallel::detectCores()) {
        warning("numCores is greater than number of cores available. Set numCores as maximum number of cores available.")
    }

    numCores <- min(numCores, parallel::detectCores())
    message(as.character(numCores), " cores are used for parallel computation.")

    if (parallelType == "MulticoreParam") {
        parallelParam <- MulticoreParam(workers = numCores)

    } else if (parallelType == "SnowParam") {
        parallelParam <- SnowParam(workers = numCores)
    } else {
        stop("'--parallelType' should be 'MulticoreParam' or 'SnowParam'.")
    }
    Params$QCMetrics$BPPARAM <- parallelParam
    Params$emptyDrops$BPPARAM <- parallelParam
    Params$doubletFinder$nCores <- numCores

}

# Helper function to remove "extension" from gene/transcript names after running kb python
# https://github.com/satijalab/seurat/issues/1049
RenameGenesSeurat <- function(obj = NULL, newnames = NULL) { # Replace gene names in obj@assays$RNA@counts, @data and @scale.data.
  RNA <- obj@assays$RNA  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

# Modify object
modifySCE <- function(seuratObj) {
  renamed <- gsub(x = rownames(seuratObj), pattern = "\\..*$", replacement = "")
  new <- RenameGenesSeurat(obj = seuratObj, newnames = renamed)
  sce <- convertSeuratToSCE(new)
  return(sce)
  
}
# read RDS object
seu <- readRDS(rds_in)
#Extract the spliced assay if necessary
if (isvelo) {
  seu <- seu$sf
}
#Set sample col
sample_col <- ifelse(replicates,"technical_replicates","descriptive_name")
#Modify sce object
sce <- modifySCE(seu)
rm(seu)
# Select QC algorithms
cellQCAlgos = c("QCMetrics", "scDblFinder")
collectionName = NULL
# Run cell QC algorithms
if (tolower(data_type) == "cell") {
    message(paste0(date(), " .. Running cell QC"))
    #Import mitochondrial gene collection
    if (isTRUE(detect_mito)) {
      #Import mitoset
      mitoset <- strsplit(mito_set,"-")
      subset_name <- stringr::str_to_title(mitoset[[1]])
      subset_name <- paste(c(subset_name, 'Mito'), collapse='')
      collectionName = subset_name
      sce <- importMitoGeneSet(sce, reference = mitoset[[1]][1], id = mitoset[[1]][2], by = "rownames", collectionName = collectionName) 
    }
    # Run QC with mitochondrial gene collection
    cellSCE <- runCellQC(sce, sample = sce[[sample_col]],
                         algorithms = cellQCAlgos,
                         collectionName = collectionName,
                         geneSetListLocation = "rownames",
                         paramsList=Params)
    # Get UMAP
    cellSCE  <- getUMAP(inSCE = cellSCE, reducedDimName = "QC_UMAP")
}
# Run droplet QC algorithms
if (tolower(data_type) == "droplet") {
    message(paste0(date(), " .. Running droplet QC"))
    dropletSCE <- runDropletQC(inSCE = sce, sample = sce[[sample_col]], paramsList=Params)
    if (isTRUE(detect_cell)) {
        if (cell_calling == "EmptyDrops") {
            ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & dropletSCE$dropletUtils_emptyDrops_fdr < 0.01
        } else if (cell_calling == "Knee") {
            ix <- dropletSCE$dropletUtils_BarcodeRank_Knee == 1
        } else {
            ix <- dropletSCE$dropletUtils_BarcodeRank_Inflection == 1
        }
        # Needs filtering of meta Data (runCellQC) sample column
        cellSCE <- dropletSCE[,ix]
        #sample_col <- cellSCE$technical_replicates
        message(paste0(date(), " .. Running cell QC"))
        #Detect mitochondrial genes
        if (isTRUE(detect_mito)) {
          #Import mitoset
          mitoset <- strsplit(mito_set,"-")
          subset_name <- stringr::str_to_title(mitoset[[1]])
          subset_name <- paste(c(subset_name, 'Mito'), collapse='')
          collectionName = subset_name
          cellSCE <- importMitoGeneSet(cellSCE, reference = mitoset[[1]][1], id = mitoset[[1]][2], by = "rownames", collectionName = collectionName)    
        }
        # Run QC with mitochondrial gene collection
       cellSCE <- runCellQC(cellSCE, sample = cellSCE[[sample_col]],
                         algorithms = cellQCAlgos,
                         collectionName = collectionName,
                         geneSetListLocation = "rownames",
                         paramsList=Params)
        #Get UMAP
        cellSCE  <- getUMAP(inSCE = cellSCE, reducedDimName = "QC_UMAP")
  }
}

#Merge result objects
if (tolower(data_type) == "cell") {
  mergedFilteredSCE <- cellSCE
  # Generate report
  message(paste0(date(), " .. Generating cell QC report"))
  reportCellQC(inSCE = mergedFilteredSCE, output_dir = out_dir, output_file = "SCTK_CellQC.html")
  #Generate QC summary
  QCsummary <- sampleSummaryStats(mergedFilteredSCE, simple=FALSE, sample = mergedFilteredSCE[[sample_col]])
  write.csv(QCsummary, qc_out)
  # Save final rds objects
  message(paste0(date(), " .. Exporting to rds format"))
  seu.processed <- as.Seurat(mergedFilteredSCE, data = NULL)
  saveRDS(seu.processed, file = rds_out)
}

#Merge Droplet sce
if (tolower(data_type) == "droplet") {
  if (isTRUE(detect_cell)) {
    mergedDropletSCE <- mergeSCEColData(dropletSCE, cellSCE)
    mergedFilteredSCE <- mergeSCEColData(cellSCE, dropletSCE)
    #Generate Report
    message(paste0(date(), " .. Generating droplet QC report"))
    pdf(pdf_out)
    print(plotEmptyDropsResults(
      inSCE = mergedDropletSCE,
      axisLabelSize = 20,
      sample = mergedDropletSCE$technical_replicates,
      fdrCutoff = 0.01,
      dotSize = 0.5,
      defaultTheme = TRUE
    ))
    # Plot barcode rank scatter
    print(plotBarcodeRankScatter(inSCE =  mergedDropletSCE,,
                                 title = "BarcodeRanks Rank Plot",
                                 legendSize = 14))
    dev.off()
    #Generate HTML report for dropletQC
    reportDropletQC(inSCE = mergedDropletSCE, output_dir = out_dir, output_file = "SCTK_DropletQC.html")
    # Generate Cell report
    message(paste0(date(), " .. Generating cell QC report"))
    reportCellQC(inSCE = mergedFilteredSCE, output_dir = out_dir, output_file = "SCTK_CellQC.html")
    #Generate QC summary
    QCsummary <- sampleSummaryStats(mergedFilteredSCE, simple=FALSE, sample = cellSCE[[sample_col]])
    write.csv(QCsummary, qc_out) 
    # Generate final rds objects
    message(paste0(date(), " .. Exporting to rds format"))
    seu_objs.processed <- c(as.Seurat(mergedFilteredSCE, data=NULL), as.Seurat(mergedDropletSCE, data=NULL))
    names(seu_objs.processed) <- c("mergedFilteredSeu","mergedDropletSeu")
    saveRDS(seu_objs.processed, file = rds_out)  
  } else {
    mergedDropletSCE <- dropletSCE
    #Generate Report
    message(paste0(date(), " .. Generating droplet QC report"))
    pdf(pdf_out)
    print(plotEmptyDropsResults(
      inSCE = mergedDropletSCE,
      axisLabelSize = 20,
      sample = mergedDropletSCE$technical_replicates,
      fdrCutoff = 0.01,
      dotSize = 0.5,
      defaultTheme = TRUE
    ))
    # Plot barcode rank scatter
    print(plotBarcodeRankScatter(inSCE =  mergedDropletSCE,,
                                 title = "BarcodeRanks Rank Plot",
                                 legendSize = 14))
    dev.off()
    reportDropletQC(inSCE = mergedDropletSCE, output_dir = out_dir, output_file = "SCTK_DropletQC.html")
    # Save final rds objects
    message(paste0(date(), " .. Exporting to rds format"))
    seu.processed <- as.Seurat(mergedFilteredSCE, data = NULL)
    saveRDS(seu.processed, file = rds_out)
  }
}













