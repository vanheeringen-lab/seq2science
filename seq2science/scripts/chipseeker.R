# log all console output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat('# variables used for this analysis:\n')
cat('narrowpeaks      <-',   snakemake@input$narrowpeaks, '\n')
cat('gtf              <-',   snakemake@params$gtf[[1]],   '\n')
cat('\n')


suppressMessages({
    library("ChIPseeker")
    library("GenomicFeatures")
})

# load the gtf file and make it a txdb
txdb_from_gtf <- makeTxDbFromGFF(snakemake@params$gtf[[1]])

peaks_list = list()
for (narrow_peak_file in snakemake@input$narrowpeaks){
    sample_name = basename(narrow_peak_file)
    peaks_list[[sample_name]] = readPeakFile(narrow_peak_file)
}

peak_anno_list <- lapply(peaks_list, annotatePeak, TxDb=txdb_from_gtf,
                         tssRegion=c(-3000, 3000), verbose=FALSE)


png(filename=snakemake@output$img1[[1]], units = 'cm', width = 20, height = 10, res = 300)
plotAnnoBar(peak_anno_list)
dev.off()

png(filename=snakemake@output$img2[[1]],  units = 'cm', width = 20, height = 10, res = 300)
plotDistToTSS(peak_anno_list)
dev.off()
