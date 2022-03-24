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
seu_min_cells <- snakemake@config$seurat_object$min_cells
seu_min_features <- snakemake@config$seurat_object$min_features

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
cat('replicates       <- "', replicates,       '"\n', sep = "")
cat('iscite           <- "', iscite,           '"\n', sep = "")
cat('isvelo           <- "', isvelo,           '"\n', sep = "")
cat('iskite           <- "', iskite,           '"\n', sep = "")

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

#Parse sample_sheet
sample_sheet <- parse_samples(samples_tsv, genome, replicates)
# Create Seurat objects based on cite input arguments and set assay
if (quantifier == "citeseqcount") {
  mat <- read_cite_output(dir=count_dir)
  meta <- prep_cell_meta(sample, sample_sheet, colnames(mat))
  # Create SingleCellExperiment object
  sce <-
    SingleCellExperiment(
      assays = list(counts = mat),
      mainExpName = sample,
      colData = meta
    )
  saveRDS(sce, file = rds)
} 
# Create Seurat objects based on input kb workflow argument and set assay
if (quantifier == "kallistobus") {
  # kb count with '--workflow kite' parameter
  if (iskite) {
    mat <- read_count_output(dir=count_dir, name="cells_x_features")
    meta <- prep_cell_meta(sample, sample_sheet, colnames(mat))
    # Create SingleCellExperiment object
    sce <-
      SingleCellExperiment(
        assays = list(counts = mat),
        mainExpName = sample,
        colData = meta
    )
    saveRDS(sce, file = rds)
    # kb count with '--workflow Lamanno'
  } else if (isvelo) {
    mat.sf <- read_count_output(count_dir, name="spliced")
    mat.uf <- read_count_output(count_dir, name="unspliced")
    meta <- prep_cell_meta(sample, sample_sheet, colnames(mat.sf))
    sce <-
      SingleCellExperiment(
        assays = list(spliced = mat.sf, unspliced = mat.uf),
        mainExpName = sample,
        colData = meta
    )
    saveRDS(sce, file = rds)
    # kb count without '--workflow' argument   
  } else {
    mat <- read_count_output(dir=count_dir, name="cells_x_genes")
    meta <- prep_cell_meta(sample, sample_sheet, colnames(mat))
    # Create SingleCellExperiment object
    sce <-
      SingleCellExperiment(
        assays = list(counts = mat),
        mainExpName = sample,
        colData = meta
    )
    saveRDS(sce, file = rds)
  }
}
