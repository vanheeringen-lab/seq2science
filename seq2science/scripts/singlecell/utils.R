suppressMessages({
    library(singleCellTK)
})

# Export a SCE object to a variety of data types
exportSCEObjs <- function(sce, out_dir, prefix = "sample") {
    # Write SCE to RDS objects
    saveRDS(object = sce, file = file.path(out_dir, "export", "R", paste0(prefix, "_", "SCE.RDS")))
    # Export to SCE to Seurat object
    exportSCEToSeurat(sce, outputDir = file.path(out_dir, "export", "R"), prefix = prefix, copyColData = TRUE, copyDecontX = TRUE, overwrite = TRUE)
    # Export SCE to FlatFile
    exportSCEtoFlatFile(sce, outputDir = file.path(out_dir, "export", "FlatFile"), prefix = prefix)
    # Generate some summary stats
    QCsummary <- sampleSummaryStats(sce, simple = FALSE, sample = NULL)
    write.csv(QCsummary, file.path(out_dir, paste0("SCE", "_", prefix, "_summary.csv")), quote = FALSE)
}

# Modifies a SingleCellExperiment object by removing the numerical transcript suffix
modifySCE <- function(sce) {
    newnames <- gsub(x = rownames(sce), pattern = "\\..*$", replacement = "")
    if (nrow(sce) == length(newnames)) {
        rownames(sce) <- newnames
    } else {
        message(paste0(date(), " .. Unequal gene sets: nrow(sce) != nrow(newnames)"))
        quit(status = 1, save = 0)
    }
    return(sce)
}

# Create scatter plot from alternative experiments
plotAltExps <- function(out_dir, sce) {
    pdf(file.path(out_dir, "SCTK_altexps.pdf", fsep = "/"))
    for (n in altExpNames(sce)) {
        x <- paste0("altexps_", n, "_percent")
        y <- "detected"
        print(plotColData(sce, x = x, y = y))
    }
    dev.off()
}
