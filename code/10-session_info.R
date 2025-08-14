x <- c(
    "dplyr",
    "tidyr",
    "limma",
    "scran",
    "scater",
    "scuttle",
    "splatter",
    "matrixStats",
    "SingleCellExperiment",
    "edgeR",
    "ggplot2",
    "ggrastr",
    "patchwork",
    "ComplexUpset",
    "ComplexHeatmap",
    "data.table",
    "viridis",
    "muscat")

# install dependencies
if (!require(BiocManager))
    install.packages("BiocManager")
for (. in x)
    if (!require(basename(.), character.only = TRUE))
        BiocManager::install(., ask = FALSE, update = TRUE)

# capture session
for (. in x) {
    . <- gsub(".*/", "", .)
    suppressPackageStartupMessages(
        library(., character.only = TRUE))
}
si <- capture.output(sessionInfo())
writeLines(si, args[[1]])