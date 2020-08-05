GTF_file_location    <- snakemake@input$gtf
Narrow_peak_files    <- snakemake@input$narrowpeaks
output_file_1 <- snakemake@output$picture_1
output_file_2 <- snakemake@output$picture_2

suppressMessages({
library("ChIPseeker")
library("GenomicFeatures")})

#Load the GTF file and make it a TxDB
TxDB_from_GTF <- makeTxDbFromGFF(GTF_file_location)

peaks_list = list()
for (narrow_peak_file in Narrow_peak_files){
    sample_name = basename(narrow_peak_file)
    peaks_list[[sample_name]] = readPeakFile(narrow_peak_file)
}

peakAnnoList <- lapply(peaks_list, annotatePeak, TxDb=TxDB_from_GTF,
                       tssRegion=c(-3000, 3000), verbose=FALSE)


png(filename=output_file_1, units = 'cm', width = 20, height = 10, res = 300)
plotAnnoBar(peakAnnoList)
dev.off()

png(filename=output_file_2,  units = 'cm', width = 20, height = 10, res = 300)
plotDistToTSS(peakAnnoList)
dev.off()