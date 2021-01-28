suppressMessages({
    library("ChIPseeker")
    library("GenomicFeatures")
})

# snakemake variables
narrowpeaks       <- snakemake@input$narrowpeaks
assembly          <- snakemake@wildcards$assembly
gtf               <- snakemake@input$gtf[[1]]
descriptive_names <- strsplit(snakemake@params$names, "\\s+")[[1]]
out1              <- snakemake@output$img1[[1]]
out2              <- snakemake@output$img2[[1]]
log_file          <- snakemake@log[[1]]

# log all console output
log <- file(log_file, open="wt")
sink(log)
sink(log, type="message")

# log all variables for debugging purposes
cat('# variables used for this analysis:\n')
cat('narrowpeaks      <-',   narrowpeaks,       '\n')
cat('assembly         <-',   assembly,          '\n')
cat('gtf              <-',   gtf,               '\n')
cat('descriptive name <-',   descriptive_names, '\n')
cat('out1             <-',   out1,              '\n')
cat('out2             <-',   out2,              '\n')
cat('\n')

cat('Sessioninfo:\n')
sessionInfo()
cat('\n')


# load the gtf file and make it a txdb
txdb_from_gtf <- makeTxDbFromGFF(gtf)

peaks_list = list()
for (i in seq_along(narrowpeaks)) {
    if (i > length(descriptive_names) || descriptive_names[[i]] == ''){
        sample_name <- narrowpeaks[[i]]
        sample_name <- gsub(".+-([^-]+)_summits\\.bed","\\1", sample_name)
    }
    else{
        sample_name = descriptive_names[[i]]
    }
    peaks_list[[sample_name]] = readPeakFile(narrowpeaks[[i]])
}


peak_anno_list <- lapply(peaks_list, annotatePeak, TxDb=txdb_from_gtf,
                         tssRegion=c(-3000, 3000), verbose=FALSE)


fig_height = 2 + 2 * length(narrowpeaks)

png(filename=out1, units = 'cm', width = 20, height = fig_height, res = 300)
plotAnnoBar(peak_anno_list, title="")
dev.off()

png(filename=out2,  units = 'cm', width = 20, height = fig_height, res = 300)
plotDistToTSS(peak_anno_list, title="")
dev.off()
