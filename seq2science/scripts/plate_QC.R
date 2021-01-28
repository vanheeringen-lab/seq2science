QC_ERCC_384plot <- function(file, name){
  library(RColorBrewer)
  library(plot3D)
  low <- "#313695";mid="#FFFFBF"; high="#A50026"   ## three-class RdYlBu, 11 levels
  RdYlBu.orig <- colorRampPalette(c(low,mid,high))(91)
  palette <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(10) # pick which palette for plate plotting
  
  coordinates<-expand.grid(seq(1,24),rev(seq(1,16)))
  #plot(expand.grid(x = c(1:24), y=c(1:16)),main=name,ylab=NA,xlab=NA, xlim=c(1,24), cex=1.5, xaxt="n", yaxt="n") #plate layout
  plot(expand.grid(x = c(1:24), y=c(1:16)),ylab=NA,xlab=NA, xlim=c(1,24), cex=1.5, xaxt="n", yaxt="n") #plate layout
  title(name, line = 2.5)
  axis(2, at=c(1:16),labels =rev(LETTERS[1:16]), las=2)
  axis(3, at=c(1:24), labels = c(1:24))
  points(coordinates,pch=19,col=palette[cut(log10(colSums(file)),include.lowest = T,
                                            breaks=unique(c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,max(colSums(file)))))], cex=1.4)
  colkey(col=RdYlBu.orig, side=1, add=T, clab="number of ERCC-reads in log10",
         line.clab = 1.25, length=0.8, clim=c(0,2))
}

QC_umis_384plot<-function(file, name){
  library(RColorBrewer)
  library(plot3D)
  low <- "#313695";mid="#FFFFBF"; high="#A50026"   ## three-class RdYlBu, 11 levels
  RdYlBu.orig <- colorRampPalette(c(low,mid,high))(91)
  palette <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(10) # pick which palette for plate plotting
  
  coordinates<-expand.grid(seq(1,24),rev(seq(1,16)))
  UMIsFig <- plot(expand.grid(x = c(1:24), y=c(1:16)),ylab=NA,xlab=NA, xlim=c(1,24), cex=1.5, xaxt="n", yaxt="n") #plate layout
  title(name, line = 2.5)
  axis(2, at=c(1:16),labels =rev(LETTERS[1:16]), las=2)
  axis(3, at=c(1:24), labels = c(1:24))
  print(paste0("Maximum column sum of plate ",name, " :" , max(colSums(file))))
  points(coordinates,pch=19,col=palette[cut(log10(colSums(file)),
                                            breaks=unique(c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,max(colSums(file)))))], cex=1.4) # plot total non-ERCC reads/cell over layout
  # cut is used here to  separate and fall the values into different internals (10 is set to tell R 10 internals are set.)
  colkey(col=RdYlBu.orig, side=1, add=T, clab="number of unique counts in log10",
         line.clab = 1.25, length=0.8, clim=c(0,5))
  
}

plate_qc <- function(data.df = NULL, spliced.data.df=NULL, barcode_file = NULL, out.file= NULL ) {
  plate_order <- read.table(barcode_file, sep = "\t", col.names = c("well","barcode"))
  #
  # # Make a vector with all plate numbers
  platenrs <- unique(gsub("([^_]*)$", "", colnames(data.df)))
  pdf(out.file, paper = "USr")
  # # settings for the plate diagnostics pdf
  par(mfrow=c(2,2), mar = c(5,4,4,2) + 0.1, cex.main = 1)
  #
  # # Iterate over all plates, order cells in the order of visualization
  for (plate in platenrs){
    #   # use the order of cells from te barcode file (this is A1, A2, A3, etc to P24)
    primer_order <- paste(plate, plate_order$well, sep="")
    #
    #   # if wells are missing on the plate, these need to be added and given a value of 0
    missing_wells <- primer_order[!primer_order %in% colnames(spliced.data.df)]
    cols_to_add <- data.frame(matrix(ncol = length(missing_wells), nrow = length(rownames(spliced.data.df))))
    colnames(cols_to_add) <- missing_wells
    cols_to_add[is.na(cols_to_add)] <- 0
    diag_plate <- cbind(spliced.data.df[,grep(plate, colnames(spliced.data.df))], cols_to_add)
    #   # phenodata contains same cellid entry + rowname as used in dataset
    cells_order <- colnames(diag_plate[,match(primer_order, colnames(diag_plate))])
    #
    #   # match dataset cells order with wells in the visualization
    tmp <- as.matrix(diag_plate[,cells_order])
    QC_umis_384plot(tmp, paste(plate, "UMI_QC", sep = "_"))
    QC_ERCC_384plot(tmp[grep("^ERCC", rownames(diag_plate)),], paste(plate, "ERCC_QC", sep = "_"))
    #
    rm(tmp)
  }
}

