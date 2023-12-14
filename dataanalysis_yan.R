# needed packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(yan),
        logcounts = log2(as.matrix(yan) + 1)
    ), 
    colData = ann
)

# preprocessing counts (gene filter plus scaling)
x <- t(counts(sce))
x <- gene_filter(x)
x <- log2(1+x)
x <- scaling(x)


absb <- get_absb(x)

