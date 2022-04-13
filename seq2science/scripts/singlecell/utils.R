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
    newnames <- stringr::str_remove(rownames(sce), pattern = "\\.\\d+")
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
    out.pdf <- paste0("SCTK_altexps_", altExpNames(sce), ".pdf")
    pdf(file.path(out_dir, out.pdf, fsep = "/"))
    for (n in altExpNames(sce)) {
        y <- paste0("altexps_", n, "_percent")
        x <- "detected"
        # Get values from colData
        x_val <- colData(sce)[[x]]
        y_val <- colData(sce)[[y]]
        # Plot
        p <- plotColData(sce, x = x, y = y) +
            ggtitle("detected features vs. percentage alt expression") +
            theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
            geom_hline(yintercept = median(y_val), col = "red", lwd = 0.5) +
            annotate("text", # Add text for mean
                x = median(x_val),
                y = max(y_val),
                label = paste("Median =", round(median(y_val), 2)),
                col = "red",
                size = 6
            )
        print(p)
    }
    dev.off()
}
