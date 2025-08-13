suppressPackageStartupMessages({
    library(muscat)
    library(BiocParallel)
    library(edgeR)
    library(SingleCellExperiment)
})

source(args[[1]])

x <- readRDS(args[[2]])
res <- fun(x)

saveRDS(res, args[[3]])