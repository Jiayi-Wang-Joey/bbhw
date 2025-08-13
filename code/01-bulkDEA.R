suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(edgeR)
    library(scuttle)
})

sce <- readRDS(args[[1]])
b <- sumCountsAcrossCells(sce, sce$sample_id)

# DEA
tab <- table(sce$sample_id, sce$group_id)
group <- tab[,1] == 0
mm <- model.matrix(~group)
dds <- calcNormFactors(DGEList(assay(b)))
dds <- dds[filterByExpr(dds, mm),]
dds <- estimateDisp(dds, mm)
fit <- glmQLFit(dds,mm)
res <- as.data.frame(topTags(glmQLFTest(fit), Inf))
res$predFC <- predFC(dds, mm)[row.names(res),2]

saveRDS(res, args[[2]])