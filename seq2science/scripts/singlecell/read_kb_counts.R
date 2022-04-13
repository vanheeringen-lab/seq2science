suppressMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

# Snakemake variables
scripts_dir <- snakemake@params$scripts_dir
count_dir <- dirname(snakemake@input$counts)
out_dir <- snakemake@params$outdir
sample <- snakemake@params$sample
replicates <- snakemake@params$replicates
genome <- snakemake@wildcards$assembly
quantifier <- snakemake@config$quantifier
iscite <- snakemake@params$iscite
isvelo <- snakemake@params$isvelo
iskite <- snakemake@params$iskite
log_file <- snakemake@log[[1]]
samples_tsv <- snakemake@config$samples
custom_assembly_suffix <- snakemake@config$custom_assembly_suffix
use_alt_expr <- snakemake@config$sc_preprocess$use_alt_expr
alt_exp_name <- snakemake@config$sc_preprocess$alt_exp_name
alt_exp_reg <- snakemake@config$sc_preprocess$alt_exp_reg
velo_assay <- snakemake@config$sc_preprocess$velo_assay

# Log all console output
log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

# Load utils library
deseq_utils <- file.path(scripts_dir, "deseq2", "utils.R")
scrna_utils <- file.path(scripts_dir, "singlecell", "utils.R")

source(scrna_utils)
source(deseq_utils)

# Log all variables for debugging purposes
cat("# variables used for this analysis:\n")
cat('quantifier       <- "', quantifier, '"\n', sep = "")
cat('sample_sheet     <- "', samples_tsv, '"\n', sep = "")
cat('scripts_dir      <- "', scripts_dir, '"\n', sep = "")
cat('log_file         <- "', log_file, '"\n', sep = "")
cat('out_dir          <- "', out_dir, '"\n', sep = "")
cat('count_dir        <- "', count_dir, '"\n', sep = "")
cat('sample           <- "', sample, '"\n', sep = "")
cat('genome           <- "', genome, '"\n', sep = "")
cat('custom_assembly_suffix  <- "', custom_assembly_suffix, '"\n', sep = "")
cat('replicates       <- "', replicates, '"\n', sep = "")
cat('iscite           <- "', iscite, '"\n', sep = "")
cat('isvelo           <- "', isvelo, '"\n', sep = "")
cat('iskite           <- "', iskite, '"\n', sep = "")
cat('velo_assay       <- "', velo_assay, '"\n', sep = "")
cat('use_alt_expr     <- "', use_alt_expr, '"\n', sep = "")
cat('alt_exp_name     <- "', alt_exp_name, '"\n', sep = "")
cat('alt_exp_reg      <- "', alt_exp_reg, '"\n', sep = "")

cat("\n")

# Prepare data frame with cell metadata
prep_cell_meta <- function(sample, sample_sheet, cell.names) {
  sample.meta <- sample_sheet[rownames(sample_sheet) %in% sample, ]
  sample.meta <- sample.meta[rep(seq_len(nrow(sample.meta)), each = length(cell.names)), ]
  rownames(sample.meta) <- cell.names
  sample.meta <- sample.meta[match(cell.names, rownames(sample.meta)), ]
  return(sample.meta)
}
# Read kb count output and return matrix
# https://pachterlab.github.io/kallistobustools/tutorials/kb_getting_started/R/kb_intro_2_R/
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/counts_unfiltered/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge.con <- file(paste0(dir, "/counts_unfiltered/", name, ".genes.txt"))
  bc.con <- file(paste0(dir, "/counts_unfiltered/", name, ".barcodes.txt"))
  colnames(m) <- readLines(bc.con)
  rownames(m) <- readLines(ge.con)
  # Close file connections
  close(ge.con)
  close(bc.con)
  return(m)
}

# Read cite-seq-count output
# https://hoohm.github.io/CITE-seq-Count/Reading-the-output/
read_cite_output <- function(dir = "", name = "umi_count") {
  matrix_dir <- paste0(dir, "/", name, "/")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names <- read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names <- read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V
  return(mat)
}

# Filter alt experiments from main experiment
filter_alt <- function(alt_feature_prefix, mat, alt = FALSE) {
  # Check if alternative experiments are present
  alt.feature <- grepl(alt_feature_prefix, rownames(mat))
  # Stop if the regex pattern cannot be founds
  if (length(alt.feature[alt.feature == TRUE]) == 0) {
    stop("Could not filter alt exp regex from gene names")
  }
  if (isTRUE(alt)) {
    return(mat[alt.feature, ])
  } else {
    return(mat[!alt.feature, ])
  }
}

# Remove custom assembly suffix from genome for filteirng
pattern <- paste0(".*", custom_assembly_suffix)
if (isTRUE(grepl(pattern, genome))) {
  genome <- substr(genome, 1, nchar(genome) - nchar(custom_assembly_suffix))
} else {
  message(paste0(date(), "No custom assembly suffix found, using default genome:", genome))
}

# Parse sample sheet
sample_sheet <- parse_samples(samples_tsv, genome, replicates)

# Create Seurat objects based on cite input arguments and set assay
alt_exp <- list()
assays <- list()

# Result matrices
mat.endo <- NULL
mat <- NULL

# Create SingleCellExperiment S4 object based on kb workflow parameter
if (quantifier == "kallistobus") {
  # kb count with '--workflow lamanno' parameter
  if (isvelo) {
    message(paste0(date(), " .. Preparing cell matrix from kallistobus (velocity) output!"))
    # Unspliced counts
    if (velo_assay == "spliced") {
      message(paste0(date(), " .. Importing spliced count matrix!"))
      mat <- read_count_output(count_dir, name = "spliced")
    } else if (velo_assay == "unspliced") {
      message(paste0(date(), " .. Importing unspliced count matrix!"))
      mat <- read_count_output(count_dir, name = "unspliced")
    } else {
      message(paste0(date(), " .. Importing spliced count matrix!"))
      mat <- read_count_output(count_dir, name = "spliced")
    }
    # Case for quantification/kite workflows
  } else {
    message(paste0(date(), " .. Preparing cell matrices from kallistobus (non-velocity) output!"))
    mat_name <- ifelse(iskite, "cells_x_features", "cells_x_genes")
    mat <- read_count_output(dir = count_dir, name = mat_name)
  }
}
# Citeseq count
if (quantifier == "citeseqcount") {
  message(paste0(date(), " .. Preparing cell matrices from citeseqcount output!"))
  mat <- read_cite_output(dir = count_dir)
}
# Check if alt experiments are present
if (use_alt_expr) {
  mat.endo <- filter_alt(alt_exp_reg, mat)
  mat.alt <- filter_alt(alt_exp_reg, mat, alt = TRUE)
  alt_exp[[alt_exp_name]] <- SummarizedExperiment(assay = list(counts = mat.alt))
} else {
  mat.endo <- mat
}
# Set assays
assays <- list(counts = mat.endo)

# Create final SingleCellExperiment S4 object and store assays
message(paste0(date(), " .. Creating final SingleCellExperiment object!"))
# Create sce object
sce <-
  SingleCellExperiment(
    assays = assays,
    mainExpName = sample,
    colData = prep_cell_meta(sample, sample_sheet, colnames(assays$counts)),
    altExps = alt_exp
  )
# Export SCE objects to various file formats
message(paste0(date(), " .. Exporting SCE object!"))
exportSCEObjs(sce, out_dir = out_dir, prefix = "raw")
