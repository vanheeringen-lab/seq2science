# scRNA-seq: plate-based CEL-Seq2 preprocessing #

**After mapping your data with Seq2Science (https://github.com/vanheeringen-lab/seq2science) or the KB-wrapper steps in `../mapping/`, these scripts can be used for further processing.**

Use the `.yaml` file to install the needed environment in Conda (using `mamba` makes this process go much quicker!).

Run in RStudio the R Markdown file with included R-Scripts, for:
- Loading the matrices from the KB-wrapper/Seq2Science pipeline.
- Running plate 'diagnostics' by generating plate overviews for ERCC and UMI counts.
- Filtering the healthy cells and expressed genes, normalize and optionally regress out confounding factors.
- Exploring the counfounding factors in the dataset.
- Running dimensionality reduction and visualizing the variability within the dataset.

**This processing RMD file expects a naming convention with "_" to generate a phenodata table from.**
For instance, if you have multiple variables per plate you would like to store in your R-object, for filtering subsets, labelling in visualization or regression, add these variables to the plate-names seperated by an `"_"` or `"-"` and these will be used to make a metadata table.

### Final product: 
A normalized and dimensionally reduced dataset that can be used directly for performing clustering, final UMAP visualizations and velocity analysis.
