resultsdir <- snakemake@output$qc_dir
output_file <- snakemake@output$pdf
kb.dir <- paste0(snakemake@config$result_dir,"/kallistobus")

#### Seurat settings ####
dual_barcodes <- snakemake@config$seurat$dual_col_barcodes 
MT_genes_file <-  snakemake@config$seurat$MT_genes_file
spike_in_qc <- snakemake@config$seurat$spike_in_qc
mt_qc <- snakemake@config$seurat$mt_qc
plate_qc <- snakemake@config$seurat$plate_qc
extract_phenotypes <- snakemake@config$seurat$extract_phenotypes
plate_variables <- snakemake@config$seurat$plate_variables
combined_variables_to_id <- snakemake@config$seurat$combined_variables_to_id
confounders_to_test <- snakemake@config$seurat$confounders_to_test
vars_to_regress <- snakemake@config$seurat$vars_to_regress
explore_violin <- snakemake@config$seurat$explore_violin
nHVG <- snakemake@config$seurat$nHVG
n_bc <- snakemake@config$seurat$n_bc
remove_bc <- snakemake@config$seurat$remove_bc 
mt_pct_max <- snakemake@config$seurat$mt_pct_max
ercc_pct_max <- snakemake@config$seurat$ercc_pct_max
pcs_max <- snakemake@config$seurat$pcs_max
amount_cells_expr <- snakemake@config$seurat$amount_cells_expr
total_counts_tresh <- snakemake@config$seurat$total_counts_tresh
total_feat_tresh <- snakemake@config$seurat$total_feat_tresh
gene_tresh <- snakemake@config$seurat$gene_tresh
pcs_for_overview <- snakemake@config$seurat$pcs_for_overview
lab_col <- snakemake@config$seurat$lab_col
umap_col <- snakemake@config$seurat$umap_col
new_col_pattern <- snakemake@config$seurat$new_col_pattern
old_col_pattern <- snakemake@config$seurat$old_col_pattern
filtering <- snakemake@config$seurat$filtering
subset_id <- snakemake@config$seurat$subset_id
#########################
rmd <- paste0(snakemake@config$rule_dir,"/../scripts/rmd/kb_seurat_pp.rmd")

#Render markdown file
rmarkdown::render(rmd, params = list(kb.dir = kb.dir, 
                                     resultsdir = resultsdir,
                                     barcode_file = dual_barcodes,
                                     MT_genes_file = MT_genes_file,
                                     add.qc.ERCC = spike_in_qc,
                                     add.qc.MT = mt_qc,
                                     run.plate_qc = plate_qc,
                                     plate_variables = plate_variables,
                                     combined_variables_to_id = combined_variables_to_id,
                                     confounders_to_test = confounders_to_test,
                                     vars_to_regress = vars_to_regress,
                                     explore_violin = explore_violin,
                                     filtering = filtering,
                                     subset_id = subset_id,
                                     nHVG = nHVG,
                                     n_bc = n_bc,
                                     remove_bc = remove_bc,
                                     mt_pct_max = mt_pct_max,
                                     ERCC_pct_max = ercc_pct_max,
                                     pcs_max = pcs_max,
                                     amount_cells_expr = amount_cells_expr,
                                     total_counts_tresh = total_counts_tresh,
                                     total_feat_tresh = total_feat_tresh,
                                     gene_tresh = gene_tresh,
                                     pcs_for_overview = pcs_for_overview,
                                     lab_col = lab_col,
                                     umap_col = umap_col,
                                     new_col_pattern = new_col_pattern,
                                     old_col_pattern = old_col_pattern,
                                     filtering = filtering), 
                                     output_file= output_file)