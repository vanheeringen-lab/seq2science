# This script is a lightweight adaption of the SCTK_runQC.R script for seq2science
# https://github.com/compbiomed/singleCellTK/blob/master/exec/SCTK_runQC.R

suppressMessages({
  library(singleCellTK)
  library("BiocParallel")
  library(scater)
})

# Snakemake variables
rds_in <- snakemake@input$rds_raw
out_dir <- snakemake@params$outdir
scripts_dir <- snakemake@params$scripts_dir
log_file <- snakemake@log[[1]]
sample <- snakemake@params$sample
isvelo <- snakemake@params$isvelo
replicates <- snakemake@params$replicates
data_type <- snakemake@config$sc_preprocess$sctk_data_type
mito_set <- snakemake@config$sc_preprocess$sctk_mito_set
detect_cell <- snakemake@config$sc_preprocess$sctk_detect_cell
detect_mito <- snakemake@config$sc_preprocess$sctk_detect_mito
cell_calling <- snakemake@config$sc_preprocess$sctk_cell_calling
use_alt_exp <- snakemake@config$sc_preprocess$use_alt_expr
pdf_out <- file.path(out_dir, "SCTK_DropletQC_figures.pdf", fsep = "/")
numCores <- snakemake@threads

# Log all console output
log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

# Load scrna utils
scrna_utils <- file.path(scripts_dir, "singlecell", "utils.R")
source(scrna_utils)


# Log all variables for debugging purposes
cat("# variables used for this analysis:\n")
cat('log_file         <- "', log_file, '"\n', sep = "")
cat('sample           <- "', sample, '"\n', sep = "")
cat('scripts_dir      <- "', scripts_dir, '"\n', sep = "")
cat('replicates       <- "', replicates, '"\n', sep = "")
cat('isvelo           <- "', isvelo, '"\n', sep = "")
cat('rds_in           <- "', rds_in, '"\n', sep = "")
cat('out_dir          <- "', out_dir, '"\n', sep = "")
cat('pdf_out          <- "', pdf_out, '"\n', sep = "")
cat('data_type        <- "', data_type, '"\n', sep = "")
cat('detect_mito      <- "', detect_mito, '"\n', sep = "")
cat('mito_set         <- "', mito_set, '"\n', sep = "")
cat('detect_cell      <- "', detect_cell, '"\n', sep = "")
cat('cell_calling     <- "', cell_calling, '"\n', sep = "")
cat('use_alt_exp      <- "', use_alt_exp, '"\n', sep = "")

# Setup parallel type
# https://github.com/compbiomed/singleCellTK/blob/master/exec/SCTK_runQC.R
parallelType <- "MulticoreParam"
Params <- list()

message(paste0(date(), " .. Setting MulticoreParam for parallel computation to:", numCores))
parallelParam <- MulticoreParam(workers = numCores)

# Set BPARAM for individual QC components
Params$QCMetrics$BPPARAM <- parallelParam
Params$emptyDrops$BPPARAM <- parallelParam
Params$doubletFinder$nCores <- numCores

# Generate QC stats for alternative experiments (if present)
if (isTRUE(use_alt_exp)) {
  Params$QCMetrics$use_altexps <- TRUE
}
# Read RDATA and modify raw sce object
sce <- readRDS(rds_in)

# Select QC algorithms
cellQCAlgos <- c("QCMetrics", "scDblFinder", "decontX")
collectionName <- NULL
# Run cell QC algorithms
if (tolower(data_type) == "cell") {
  message(paste0(date(), " .. Running CellQC"))
  # Import mitochondrial gene collection
  if (isTRUE(detect_mito)) {
    # Import mitoset
    mitoset <- strsplit(mito_set, "-")
    subset_name <- stringr::str_to_title(mitoset[[1]])
    subset_name <- paste(c(subset_name, "Mito"), collapse = "")
    collectionName <- subset_name
    sce <- importMitoGeneSet(sce, reference = mitoset[[1]][1], id = mitoset[[1]][2], by = "rownames", collectionName = collectionName)
  }
  # Run QC with mitochondrial gene collection
  cellSCE <- runCellQC(sce,
    sample = NULL,
    algorithms = cellQCAlgos,
    collectionName = collectionName,
    geneSetListLocation = "rownames",
    paramsList = Params
  )
}
# Run droplet QC algorithms
if (tolower(data_type) == "droplet") {
  message(paste0(date(), " .. Running DropletQC"))
  dropletSCE <- runDropletQC(inSCE = sce, sample = NULL, paramsList = Params)
  if (isTRUE(detect_cell)) {
    if (cell_calling == "EmptyDrops") {
      ix <- !is.na(dropletSCE$dropletUtils_emptyDrops_fdr) & dropletSCE$dropletUtils_emptyDrops_fdr < 0.01
    } else if (cell_calling == "Knee") {
      ix <- dropletSCE$dropletUtils_BarcodeRank_Knee == 1
    } else {
      ix <- dropletSCE$dropletUtils_BarcodeRank_Inflection == 1
    }
    # Needs filtering of meta Data (runCellQC) sample column
    cellSCE <- dropletSCE[, ix]
    message(paste0(date(), " .. Running CellQC"))
    # Detect mitochondrial genes
    if (isTRUE(detect_mito)) {
      # Import mitoset
      mitoset <- strsplit(mito_set, "-")
      subset_name <- stringr::str_to_title(mitoset[[1]])
      subset_name <- paste(c(subset_name, "Mito"), collapse = "")
      collectionName <- subset_name
      cellSCE <- importMitoGeneSet(cellSCE, reference = mitoset[[1]][1], id = mitoset[[1]][2], by = "rownames", collectionName = collectionName)
    }
    # Run QC with mitochondrial gene collection
    cellSCE <- runCellQC(cellSCE,
      sample = NULL,
      algorithms = cellQCAlgos,
      collectionName = collectionName,
      geneSetListLocation = "rownames",
      paramsList = Params
    )
  }
}

# Merge result objects
if (tolower(data_type) == "cell") {
  mergedFilteredSCE <- cellSCE
  # Generate CellQC report
  message(paste0(date(), " .. Generating CellQC report"))
  reportCellQC(inSCE = mergedFilteredSCE, output_dir = out_dir, output_file = "SCTK_CellQC.html")
  # Save final rds objects
  if (isTRUE(use_alt_exp)) {
    plotAltExps(out_dir, mergedFilteredSCE)
  }
  message(paste0(date(), " .. Exporting SCE object!"))
  exportSCEObjs(mergedFilteredSCE, out_dir = out_dir, prefix = "sctk")
}

# Merge Droplet SingleCellExperiment object
if (tolower(data_type) == "droplet") {
  if (isTRUE(detect_cell)) {
    mergedDropletSCE <- mergeSCEColData(dropletSCE, cellSCE)
    mergedFilteredSCE <- mergeSCEColData(cellSCE, dropletSCE)
    # Generate DropletQC Report
    message(paste0(date(), " .. Generating DropletQC report"))
    pdf(file.path(out_dir, "SCTK_DropletQC_figures.pdf", fsep = "/"))
    print(plotEmptyDropsResults(
      inSCE = mergedDropletSCE,
      axisLabelSize = 20,
      sample = NULL,
      fdrCutoff = 0.01,
      dotSize = 0.5,
      defaultTheme = TRUE
    ))
    # Plot barcode rank scatter
    print(plotBarcodeRankScatter(
      inSCE = mergedDropletSCE, ,
      title = "BarcodeRanks Rank Plot",
      legendSize = 14
    ))
    dev.off()
    # Generate HTML report for DropletQC
    reportDropletQC(inSCE = mergedDropletSCE, output_dir = out_dir, output_file = "SCTK_DropletQC.html")
    # Generate CellQC report
    message(paste0(date(), " .. Generating CellQC report"))
    reportCellQC(inSCE = mergedFilteredSCE, output_dir = out_dir, output_file = "SCTK_CellQC.html")
    # Generate report for alternative experiments
    if (isTRUE(use_alt_exp)) {
      plotAltExps(out_dir, mergedFilteredSCE)
    }
    # Generate final RDATA object
    message(paste0(date(), " .. Exporting SCE object!"))
    # sce.processed <- list(cellsce = mergedFilteredSCE, dropletsce = mergedDropletSCE)
    exportSCEObjs(mergedFilteredSCE, out_dir = out_dir, prefix = "sctk")
  } else {
    mergedDropletSCE <- dropletSCE
    # Generate DropletQC Report
    message(paste0(date(), " .. Generating DropletQC report"))
    pdf(file.path(out_dir, "SCTK_DropletQC_figures.pdf", fsep = "/"))
    print(plotEmptyDropsResults(
      inSCE = mergedDropletSCE,
      axisLabelSize = 20,
      sample = NULL,
      fdrCutoff = 0.01,
      dotSize = 0.5,
      defaultTheme = TRUE
    ))
    # Plot barcode rank scatter
    print(plotBarcodeRankScatter(
      inSCE = mergedDropletSCE,
      title = "BarcodeRanks Rank Plot",
      legendSize = 14
    ))
    dev.off()
    # Create DropletQC report
    reportDropletQC(inSCE = mergedDropletSCE, output_dir = out_dir, output_file = "SCTK_DropletQC.html")
    # Export SCE objects
    message(paste0(date(), " .. Exporting SCE object!"))
    exportSCEObjs(mergedDropletSCE, out_dir = out_dir, prefix = "sctk")
  }
}
