suppressPackageStartupMessages({
    library(muscat)
    library(BiocParallel)
    library(edgeR)
    library(SingleCellExperiment)
})

fun <- \(x) {
    sz <- c(2,3,4,5,6,7,8)
    pb <- aggregateData(x, by = c("cluster_id", "sample_id"))
    res <- bplapply(sz, 
                      BPPARAM=MulticoreParam(6, RNGseed=1234, progress=FALSE), 
                      \(s){
                    
        sw <- c(sample(which(pb$group_id=="A"),s),sample(which(pb$group_id=="B"),s))
        pb <- pb[,sw]
        d <- lapply(setNames(assayNames(pb),assayNames(pb)), FUN=\(x){
            tryCatch({
                assays(pb) <- assays(pb)[x]
                mm <- model.matrix(~group_id, data=as.data.frame(colData(pb)))
                dds <- calcNormFactors(DGEList(assay(pb)))
                dds <- dds[filterByExpr(dds,mm),]
                dds <- estimateDisp(dds, mm)
                fit <- glmQLFit(dds,mm)
                d <- as.data.frame(topTags(glmQLFTest(fit), Inf))
                d$gene <- row.names(d)
                d$pbDEA_size <- paste0(s,"v", s)
                d
            }, error=\(e){ warning(e); NULL })
        })
        dplyr::bind_rows(d[!sapply(d, is.null)], .id="celltype")
    })
    names(res) <- assayNames(pb)
    res
}
