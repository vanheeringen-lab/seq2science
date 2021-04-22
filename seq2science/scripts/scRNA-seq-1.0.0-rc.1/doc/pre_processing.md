# scRNA-seq # 
Preprocessing datasets, starting from mapping the raw data from the plate-based CEL-seq2 protocol (`mapping/`) and generating a processed and QC-ed dataset for downstream applications (`processing/`). 

The dataset will be preprocessed for use of RNA Velocity, therefore also processing the unspliced reads.

- From fastq to preprocessed counttable (for in-house CELSeq2 method), with Kallisto | Bustools wrapper.
- QC-ed, normalized and dimensionality reduced in custom R scripts, with Scater and Seurat. 

In the `genome_addons/`, various files are stored useful for genome building with KB-wrapper (spike-ins and reporter protein sequences).
