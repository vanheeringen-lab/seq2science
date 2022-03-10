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
data_type <- snakemake@config$sctk$data_type
mito_set <- snakemake@config$sctk$mito_set
rds_out <- file.path(out_dir, "seu_obj_processed.RData",    fsep="/" )
qc_out <-  file.path(out_dir, "SCTK_cellQC_summary.csv",    fsep="/" )
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
cat('rds_in           <- "', rds_in,           '"\n', sep = "")
cat('out_dir          <- "', out_dir,          '"\n', sep = "")
cat('rds_out          <- "', out_dir,          '"\n', sep = "")
cat('qc_out           <- "', qc_out,           '"\n', sep = "")
cat('pdf_out          <- "', pdf_out,          '"\n', sep = "")
cat('data_type        <- "', data_type,        '"\n', sep = "")
cat('mito_set          <- "', mito_set,          '"\n', sep = "")


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
    #Params$doubletCells$BPPARAM <- parallelParam
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
sce <- modifySCE(seu)
#Perform filtering steps for Droplet based methods
if(tolower(data_type) == "droplet") {
  sce <-
    runDropletQC(
      sce,
      algorithms = c("QCMetrics", "emptyDrops", "barcodeRanks"),
      sample = sce$technical_replicates,
      paramsList = Params
    )
  reportDropletQC(inSCE = sce, output_dir = out_dir, output_file = "DropletQC.html")
  # Filtering
  sce <-
    subsetSCECols(sce, colData = 'dropletUtils_BarcodeRank_Inflection == 1')
  sce <-
    subsetSCECols(sce, colData = '!is.na(sce$dropletUtils_emptyDrops_fdr)')
  sce <-
    subsetSCECols(sce, colData = 'sce$dropletUtils_emptyDrops_fdr < 0.01')
  # Plot Results
  pdf(pdf_out)
  print(plotEmptyDropsResults(
    inSCE = sce,
    axisLabelSize = 20,
    sample = NULL,
    fdrCutoff = 0.01,
    dotSize = 0.5,
    defaultTheme = TRUE
  ))
  # Plot barcode rank scatter
  print(plotBarcodeRankScatter(inSCE = sce,
                         title = "BarcodeRanks Rank Plot",
                         legendSize = 14))
  #Generate HTML report for droplet
  dev.off()
}

#Import Mito Gene set for Quality contro
mitoset <- strsplit(mito_set,"-")[[1]]
sce <- importMitoGeneSet(sce, reference = mitoset[1], id = mitoset[2], by = "rownames", collectionName = "mito")
sce <- runCellQC(sce, sample = sce$technical_replicates,
                   algorithms = c("QCMetrics", "scDblFinder"),
                   collectionName = "mito",
                   geneSetListLocation = "rownames",
                   paramsList=Params)
sce <- getUMAP(inSCE = sce, reducedDimName = "QC_UMAP")
#Generate HTML report for Cell analysis
reportCellQC(sce, output_dir = out_dir, output_file = "CellQC.html")

# Generate summary
QCsummary <- sampleSummaryStats(sce, simple=FALSE, sample = sce$technical_replicates)
write.csv(QCsummary, qc_out)

#Save final Seurat object
seu.processed <- as.Seurat(sce, data = NULL)
saveRDS(seu.processed, file = rds_out)


