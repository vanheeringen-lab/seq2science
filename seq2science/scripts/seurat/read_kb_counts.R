suppressMessages({
    library(Matrix)
    library(Seurat)
})

# Snakemake variables
count_dir <- dirname(snakemake@input$counts)
rds <- snakemake@output[[1]] 
sample <- snakemake@params$sample
replicates <- snakemake@params$replicates
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

# Read samplesheet
sample_sheet <- read.table(samples_tsv, sep = '\t', header = TRUE)

# Log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('sample_sheet <- "', samples_tsv,     '"\n', sep = "")
cat('log_file     <- "', log_file,        '"\n', sep = "")
cat('count_dir       <- "', count_dir,          '"\n', sep = "")
cat('rds          <- "', rds,             '"\n', sep = "")
cat('sample       <- "', sample,          '"\n', sep = "")
cat('replicates       <- "', replicats,          '"\n', sep = "")
cat('iscite       <- "', iscite,          '"\n', sep = "")
cat('isvelo       <- "', isvelo,          '"\n', sep = "")
cat('iskite       <- "', iskite,          '"\n', sep = "")
cat('seu_min_cells  <- "', seu_min_cells, '"\n', sep = "")
cat('seu_min_features <- "', seu_min_features,  '"\n', sep = "")


cat('\n')

# Read cell meta data from samplesheet
prep_cell_meta <- function(seu, sample_sheet) {
  blacklist <- c("descriptive_name", "technical_replicates")
  meta <- colnames(sample_sheet[setdiff(names(sample_sheet),blacklist)])
  sample.meta <- sample_sheet[sample_sheet$sample==seu@project.name,]
  sample.meta <- sample.meta[meta]
  sample.meta <- sample.meta[rep(seq_len(nrow(sample.meta)), each = ncol(seu)),]
  rownames(sample.meta) <- rownames(FetchData(seu,"ident"))
  sample.meta <- sample.meta[match(rownames(FetchData(seu,"ident")),rownames(sample.meta)),]
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

# Create Seurat objects based on cite input arguments and set assay
if (quantifier == "citeseqcount") {
  seu <- CreateSeuratObject(counts = read_cite_output(dir=count_dir), assay = "ADT", project = sample, 
                                                      min.cells = seu_min_cells, min.features = seu_min_features)
  if (!replicates) {
    seu@meta.data <- prep_cell_meta(seu, sample_sheet)
  }
  saveRDS(seu, file = rds)  
} 
# Create Seurat objects based on input kb workflow argument and set assay
if (quantifier == "kallistobus") {
  # kb count with '--workflow kite' parameter
  if (iskite) {
    seu <- CreateSeuratObject(counts = read_count_output(count_dir, name="cells_x_features"), assay = "ADT", project = sample, 
                                                       min.cells = seu_min_cells, min.features = seu_min_features)
    if (!replicates)) {
      seu@meta.data <- prep_cell_meta(seu, sample_sheet)
    }
    saveRDS(seu, file = rds)
    # kb count with '--workflow Lamanno'
  } else if (isvelo) {
    seu.sf <- CreateSeuratObject(counts = read_count_output(count_dir, name="spliced"), assay = "sf", project = sample, 
                                                            min.cells = seu_min_cells, min.features = seu_min_features)
    seu.uf <- CreateSeuratObject(counts = read_count_output(count_dir, name="unspliced"), assay = "uf", project = sample, 
                                                            min.cells = seu_min_cells, min.features = seu_min_features)
    if (!replicates) {
      seu.sf@meta.data <- prep_cell_meta(seu.sf, sample_sheet)
      seu.uf@meta.data <- prep_cell_meta(seu.uf, sample_sheet)
    }
    seu_objs <- c(seu.sf, seu.uf)
    names(seu_objs) <- c("sf","uf")
    saveRDS(seu_objs, file = rds) 
    # kb count without '--workflow' argument   
  } else {
    seu <- CreateSeuratObject(counts = read_count_output(count_dir, name="cells_x_genes"), assay = "RNA", project = sample, 
                                                         min.cells = seu_min_cells, min.features = seu_min_features)
    
    if (replicates)) {
      seu@meta.data <- prep_cell_meta(seu, sample_sheet)
    }
    saveRDS(seu, file = rds)
  }
}
