#' parse the design contrast
#' a contrast is always in the form 'batch+condition_group1_group2', where batch(+) is optional
parse_contrast <- function(contrast) {
  batch <- NA
  contr <- contrast

  # batch name
  if (grepl('\\+', contrast)) {
    batch <- strsplit(contrast, '\\+')[[1]][1]
    contr <- strsplit(contrast, '\\+')[[1]][2]
  }

  # group names
  groups <- strsplit(contr, '_')[[1]]
  groups <- tail(groups,2)

  # column name
  n <- gregexpr(pattern=paste0("_", groups[1], "_", groups[2]), contr)[[1]][1] -1
  condition <- substr(contr, 1, n)

  return(
    list(
      "batch"=batch,
      "condition"=condition,
      "groups"=groups
    )
  )
}


#' load the samples.tsv file into a dataframe,
#' filter for the right assembly,
#' collapse technical replicates (if present) and
#' convert to descriptive names (if present).
parse_samples <- function(samples_file, assembly_name, replicates) {
  samples <- read.delim(samples_file, sep = "\t", na.strings = "", comment.char = "#", stringsAsFactors = F, row.names = "sample", check.names = F, colClasses="character")
  colnames(samples) <- gsub("\\s+", "", colnames(samples))  # strip whitespace from column names

  # custom assemblies have a suffix not present in the samples.tsv
  if (!(assembly_name %in% samples$assembly)){
    for (assembly in unique(samples$assembly)){
      assembly <- trimws(assembly)
      if (startsWith(assembly_name, assembly)){
        assembly_name <- assembly
        break
      }
    }
  }

  # drop all samples of other assemblies
  samples <- subset(samples, assembly == assembly_name)

  # collapse technical replicates
  # (and use these names for the counts table later)
  if ("technical_replicates" %in% colnames(samples) & isTRUE(replicates)) {
    to_rename <- is.na(samples$technical_replicates)
    samples$technical_replicates[to_rename] <- as.character(rownames(samples)[to_rename])
    samples <- subset(samples, !duplicated(technical_replicates))
    rownames(samples) <- samples$technical_replicates
  }

  # fill in blank descriptive names
  # (if no replicates are found, we can use descriptive names for the counts table later)
  if ("descriptive_name" %in% colnames(samples)) {
    to_rename <- is.na(samples$descriptive_name)
    samples$descriptive_name[to_rename] <- as.character(rownames(samples)[to_rename])
    rownames(samples) <- samples$descriptive_name
  }

  if (nrow(samples) == 0){
    stop("Something went wrong filtering the samples file! No samples remaining.")
  }

  return(samples)
}


#' run DESeq2
run_deseq2 <- function(counts, coldata, design, threads=1, single_cell=FALSE) {
  parallel <- FALSE
  if (threads > 1) {
    BiocParallel::register(
      BiocParallel::MulticoreParam(threads)
    )
    parallel <- TRUE
  }

  cat('Constructing DESeq object... \n')
  cat('Tip: errors directly below this line are most likely DESeq2 related.\n\n')
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = design,
  )
  cat('\nFinished constructing DESeq object.\n\n')

  if (single_cell) {
    # source: DESeq2
    # https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/
    # DESeq2.html#recommendations-for-single-cell-analysis

    # Model zero component using zinbwave
    isZeroInfated = FALSE  # TODO: expose? It's unclear if zero inflation is a thing
    if (isZeroInfated) {
      # low count filter - at least 10 cells with a read count of 5 or more
      keep <- rowSums(counts >= 5) >= 10
      if (sum(keep) > 10) {zinb <- dds[keep,]} else {zinb <- dds}

      # epsilon setting as recommended by the ZINB-WaVE integration paper
      zinb <- zinbwave::zinbwave(
        zinb,
        K=0,
        observationalWeights=TRUE,
        BPPARAM=BiocParallel::SerialParam(),
        epsilon=1e12,
      )

      dds <- DESeq2::DESeqDataSet(zinb, design=design)
    }

    # Estimate size factors
    scr <- scran::computeSumFactors(dds)
    DESeq2::sizeFactors(dds) <- DESeq2::sizeFactors(scr)

    # Estimate dispersion and DE
    dds <- DESeq2::DESeq(
      dds,
      test="LRT",
      reduced=~1,  # required for LRT
      useT=TRUE,
      minmu=1e-6,
      minReplicatesForReplace=Inf,
      parallel=parallel,
#       fitType = "glmGamPoi",  # silently deactivated with parallel=TRUE
    )
  } else {
    dds <- DESeq2::DESeq(dds, parallel=parallel)
  }
  cat('\n')

  return(dds)
}


#' save a diffexp table with all genes, not just the expressed genes
save_complete_diffexp <- function(resLFC, counts, diffexp_file) {
  expressed_genes <- as.data.frame(resLFC[order(resLFC$padj),])

  missing_genes <- rownames(counts)[!(rownames(counts) %in% rownames(expressed_genes))]
  if (length(missing_genes) > 0){
      unexpressed_genes <- as.data.frame(matrix(data = NA, ncol = ncol(expressed_genes), nrow = length(missing_genes)))
      rownames(unexpressed_genes) <- missing_genes
      colnames(unexpressed_genes) <- colnames(expressed_genes)
      unexpressed_genes[,'baseMean'] <- 0
      all_genes <- rbind(expressed_genes, unexpressed_genes)
  } else {
      all_genes <- expressed_genes
  }
  write.table(all_genes, file=diffexp_file, quote = F, sep = '\t', col.names=NA)
}


#' model the effect of batch correction in the correlation heatap.
#' DESeq2 applies this internally as well.
batch_corrected_vst <- function(dds) {
  nonblind_vst <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  mat <- SummarizedExperiment::assay(nonblind_vst)
  mat <- limma::removeBatchEffect(mat, nonblind_vst$batch)

  batchcorr_vst <- nonblind_vst
  SummarizedExperiment::assay(batchcorr_vst) <- mat
  batchcorr_vst <- DESeq2::DESeqTransform(batchcorr_vst)

  return(batchcorr_vst)
}


#' accepts a large DESeqDataSet, normalizes and removes the batch effects,
#' and returns a batch corrected count matrix
batch_corrected_counts <- function(dds) {
  nonblind_vst <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  mat <- SummarizedExperiment::assay(nonblind_vst)
  mat <- limma::removeBatchEffect(mat, nonblind_vst$batch)

  # mat contains the normalized and batch corrected variant stabilized data (mat = vst.fn(batch_corr_counts))
  # next, we invert this function (from DESeq2::getVarianceStabilizedData) to get normalized and batch corrected counts (batch_corr_counts)

  # vst.fn <- function(q) {
  #     log((1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] *
  #         q + 2 * sqrt(coefs["asymptDisp"] * q * (1 + coefs["extraPois"] +
  #         coefs["asymptDisp"] * q)))/(4 * coefs["asymptDisp"]))/log(2)
  # }
  coefs <- attr(DESeq2::dispersionFunction(dds), "coefficients")
  x <- mat * log(2)
  x <- exp(1)^x * (4 * coefs["asymptDisp"])
  x <- (x - (1 + coefs["extraPois"]))/2
  batch_corr_counts <- x^2 / (coefs["asymptDisp"] * (coefs["extraPois"] + 2 * x + 1))

  return(batch_corr_counts)
}


#' snippet from https://support.bioconductor.org/p/91218/ to convert counts to TPM
counts2tpm <- function(mat.counts, mat.gene_lengths) {
  x <- mat.counts / mat.gene_lengths
  mat.tpm <- t( t(x) * 1e6 / colSums(x) )
  return(mat.tpm)
}


#' heatmap aesthetics variable to the number of samples
heatmap_aesthetics <- function(num_samples){
  if (num_samples < 16) {
    cell_dimensions <- as.integer(160/num_samples)  # minimal size to fit the legend
    fontsize        <- 8
    fontsize_number <- 8
  } else if (num_samples < 24) {
    cell_dimensions <- 10  # pleasant size
    fontsize        <- 8
    fontsize_number <- 8
  } else if (num_samples < 32) {
    cell_dimensions <- as.integer(25 - 0.625*num_samples)  # linear shrink
    fontsize        <- 8 - 0.15*num_samples
    fontsize_number <- 7 - 0.15*num_samples
  } else {
    cell_dimensions <- 5  # minimal size
    fontsize        <- 3.2
    fontsize_number <- 2.5
  }
  show_colnames <- ifelse(num_samples > 28, TRUE, FALSE)
  show_rownames <- !show_colnames
  colors <- colorRampPalette( RColorBrewer::brewer.pal(9, "RdYlBu") )(255)

  return(
    list(
      "cell_dimensions"=cell_dimensions,
      "fontsize"=fontsize,
      "fontsize_number"=fontsize_number,
      "show_colnames"=show_colnames,
      "show_rownames"=show_rownames,
      "colors"=colors
    )
  )
}


#' assign names from coldata to the rows and columns of a matrix
#' uses descriptive names if available, else rownames (can be technical replicates/sample names)
heatmap_names <- function(mat, coldata) {
  names <- rownames(coldata)
  rownames(mat) <- names
  colnames(mat) <- names
  return(mat)
}


heatmap_plot <- function(mat, title, heatmap_aes, legend_aes, out_pdf) {
  #' rotate the dendogram to best match the order in the samples.tsv
  callback <- function(hc, mat){
    hc <- dendextend::rotate(hc, order=colnames(mat))
    return(hc)
  }

  pheatmap::pheatmap(
    mat,
    clustering_callback = callback,
    main = title,
    angle_col = 45,
    show_colnames = heatmap_aes$show_colnames,  # show names underneath if the image gets to wide
    show_rownames = heatmap_aes$show_rownames,
    fontsize = heatmap_aes$fontsize,
    legend_breaks = legend_aes$breaks,
    legend_labels = legend_aes$labels,
    # display_numbers = T,  # show values in the plot
    fontsize_number = heatmap_aes$fontsize_number,
    cellwidth  = heatmap_aes$cell_dimensions,
    cellheight = heatmap_aes$cell_dimensions,
    color = heatmap_aes$colors,
    filename = out_pdf
  )
}

#' create a PNG file from a PDF file
pdf2png <- function(pdf, png) {
  pdftools::pdf_convert(
    pdf,
    format = "png",
    filenames = png,
    dpi = 300,  # standard minimum
    antialias = TRUE,
    verbose = TRUE
  )
}
