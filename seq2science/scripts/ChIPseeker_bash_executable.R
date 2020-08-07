args = commandArgs(trailingOnly = TRUE)
if (length(args)<2) {
    print('Need both a GTF or GFF3 file, and minimun one narrow peak file to annotate')}
if (length(args)>1) {
    print('using the following files:')
    print((paste("supplied GTF/GFF3 File: ", args[1])))
    print((paste("using narrowPeak file:", args[2:length(args)])))}

suppressMessages({
library("ChIPseeker")
library("GenomicFeatures")})

#Load the GTF file and make it a TxDB
TxDB_from_GTF <- makeTxDbFromGFF(args[1])

peaks_list = list()
for (narrow_peak_file in args[2:length(args)]){
    print(narrow_peak_file)
    sample_name = basename(narrow_peak_file)
    peaks_list[[sample_name]] = readPeakFile(narrow_peak_file)
}

peakAnnoList <- lapply(peaks_list, annotatePeak, TxDb=TxDB_from_GTF,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

png(filename="/home/jsmits/peaks_annotated_2_TSS_mqc.png", units = 'cm', width = 20, height = 10, res = 300)

plotAnnoBar(peakAnnoList)
dev.off()

png(filename="/home/jsmits/peaks_dist_2_TSS_mqc.png", units = 'cm', width = 20, height = 10, res = 300)

plotDistToTSS(peakAnnoList)
dev.off()