# log all console output
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat('# variables used for this analysis:\n')
cat('narrowpeaks      <-',   snakemake@input$narrowpeaks, '\n')
cat('gtf              <-',   snakemake@params$gtf[[1]],   '\n')
cat('descriptive name <-',   snakemake@params$names,      '\n')
cat('\n')


suppressMessages({
    library("ChIPseeker")
    library("GenomicFeatures")
})

descriptive_names <- strsplit(snakemake@params$names, "\\s+")[[1]]

# load the gtf file and make it a txdb
txdb_from_gtf <- makeTxDbFromGFF(snakemake@params$gtf[[1]])

peaks_list = list()
for (i in seq_along(snakemake@input$narrowpeaks)) {
    sample_name = descriptive_names[[i]]
    peaks_list[[sample_name]] = readPeakFile(snakemake@input$narrowpeaks[[i]])
}

peak_anno_list <- lapply(peaks_list, annotatePeak, TxDb=txdb_from_gtf,
                         tssRegion=c(-3000, 3000), verbose=FALSE)


fig_height = 2 + 2 * length(snakemake@input$narrowpeaks)

png(filename=snakemake@output$img1[[1]], units = 'cm', width = 20, height = fig_height, res = 300)
plotAnnoBar(peak_anno_list, title="")
dev.off()

png(filename=snakemake@output$img2[[1]],  units = 'cm', width = 20, height = fig_height, res = 300)
plotDistToTSS(peak_anno_list, title="")
dev.off()
