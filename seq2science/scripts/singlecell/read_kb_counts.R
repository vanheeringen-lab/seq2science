suppressMessages({
    library(Matrix)
    library(SingleCellExperiment)
})

# Snakemake variables
scripts_dir  <- snakemake@params$scripts_dir
count_dir <- dirname(snakemake@input$counts)
rds <- snakemake@output[[1]] 
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

# Log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

#Load utils library
deseq_utils <- file.path(scripts_dir,"utils.R")
source(deseq_utils)

# Log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('quantifier       <- "', quantifier,       '"\n', sep = "")
cat('sample_sheet     <- "', samples_tsv,      '"\n', sep = "")
cat('scripts_dir      <- "', scripts_dir,      '"\n', sep = "")
cat('log_file         <- "', log_file,         '"\n', sep = "")
cat('count_dir        <- "', count_dir,        '"\n', sep = "")
cat('rds              <- "', rds,              '"\n', sep = "")
cat('sample           <- "', sample,           '"\n', sep = "")
cat('genome           <- "', genome,           '"\n', sep = "")
cat('custom_assembly_suffix  <- "', custom_assembly_suffix, '"\n', sep = "")
cat('replicates       <- "', replicates,       '"\n', sep = "")
cat('iscite           <- "', iscite,           '"\n', sep = "")
cat('isvelo           <- "', isvelo,           '"\n', sep = "")
cat('iskite           <- "', iskite,           '"\n', sep = "")
cat('use_alt_expr     <- "', use_alt_expr,     '"\n', sep = "")
cat('alt_exp_name     <- "', alt_exp_name,     '"\n', sep = "")
cat('alt_exp_reg      <- "', alt_exp_reg,      '"\n', sep = "")

cat('\n')

#Prep cell metadata
prep_cell_meta <- function(sample, sample_sheet, cell.names) {
  sample.meta <- sample_sheet[rownames(sample_sheet) %in% sample,]
  sample.meta <- sample.meta[rep(seq_len(nrow(sample.meta)), each = length(cell.names)),]
  rownames(sample.meta) <- cell.names
  sample.meta <- sample.meta[match(cell.names,rownames(sample.meta)),]
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
read_cite_output <- function(dir="", name="umi_count") {
  matrix_dir=paste0(dir,"/",name,"/")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V
  return(mat)
}

# Filter alt experiments from main experiment
filter_alt <- function(alt_feature_prefix, mat, alt = FALSE) {
  # Check if alternative experiments are present
  alt.feature <- grepl(alt_feature_prefix, rownames(mat))
  # Stop if the regex pattern cannot be founds
  if (length(alt.feature[alt.feature == TRUE]) == 0 ) {
    stop("Could not filter alt exp regex from gene names")
  }
  if (isTRUE(alt)) {
    return(mat[alt.feature,])
  }
  else {
    return(mat[!alt.feature,])
  }
}

# Remove custom assembly suffix from genome for filteirng
if (custom_assembly_suffix != "") {
  genome <- gsub(custom_assembly_suffix, "", genome)
}
# Parse sample sheet
sample_sheet <- parse_samples(samples_tsv, genome, replicates)
# Create Seurat objects based on cite input arguments and set assay
alt_exp <- list()
assays <- list()
# Create Seurat objects based on input kb workflow argument and set assay
if (quantifier == "kallistobus") {
  # kb count with '--workflow kite' parameter
  if (isvelo) {
    message(paste0(date(), " .. Preparing cell matrices from kallistobus (velocity) output!"))
    # Unspliced counts
    mat.sf.endo <- NULL
    mat.uf.endo <- NULL
    mat.sf <- read_count_output(count_dir, name="spliced")
    mat.uf <- read_count_output(count_dir, name="unspliced")
    if (use_alt_expr) {
      mat.sf.endo <- filter_alt(alt_exp_reg, mat.sf)
      mat.uf.endo <- filter_alt(alt_exp_reg, mat.uf)
      # Filter alt experiment from main experiment
      mat.sf.alt <- filter_alt(alt_exp_reg, mat.sf, alt = TRUE)
      mat.uf.alt <- filter_alt(alt_exp_reg, mat.uf, alt = TRUE)
      # Store experiments for spliced and unsplaced assays
      alt_exp[[paste0(alt_exp_name,"_sf")]] <- SummarizedExperiment(assay = list(counts=mat.sf.alt)) 
      alt_exp[[paste0(alt_exp_name,"_uf")]] <- SummarizedExperiment(assay = list(counts=mat.uf.alt)) 
    } else {
      mat.sf.endo <- mat.sf
      mat.uf.endo <- mat.uf
    }
    # create assays
    assays <- list(counts = mat.sf.endo, unspliced = mat.uf.endo)
  
  # Case for quantification/kite workflows
  } else {
    message(paste0(date(), " .. Preparing cell matrices from kallistobus (non-velocity) output!"))
    mat_name <- ifelse(iskite, "cell_x_features", "cell_x_genes")
    mat <- read_count_output(dir=count_dir, name=mat_name)
    mat.endo <- NULL
    if (use_alt_expr) {
      mat.endo <- filter_alt(alt_exp_reg, mat)
      mat.alt <- filter_alt(alt_exp_reg, mat, alt = TRUE)
      alt_exp[[alt_exp_name]] <- SummarizedExperiment(assay = list(counts=mat.alt)) 
    } else {
      mat.endo <- mat
    }
    assays <- list(counts = mat.endo)
  }    
}

if (quantifier == "citeseqcount") {
  message(paste0(date(), " .. Preparing cell matrices from citeseqcount output!"))
  mat <- read_cite_output(dir=count_dir)
  mat.endo <- NULL
  if (use_alt_expr) {
    mat.endo <- filter_alt(alt_exp_reg, mat)
    mat.alt <- filter_alt(alt_exp_reg, mat, alt = TRUE)
    alt_exp[[alt_exp_name]] <- SummarizedExperiment(assay = mat.alt) 
  } else {
    mat.endo <- mat
  }
  assays <- list(counts = mat.endo)
}
  
# Create final sce object and store main and alternative experiments (if present) 
message(paste0(date(), " .. Creating final SingleCellExperiment object!"))
sce <-
  SingleCellExperiment(
    assays = assays,
    mainExpName = sample,
    colData = prep_cell_meta(sample, sample_sheet, colnames(assays$counts)),
    altExps =  alt_exp)

# Save RDS object
message(paste0(date(), " .. Saving object to RDSD file!"))
saveRDS(sce, file = rds)  


