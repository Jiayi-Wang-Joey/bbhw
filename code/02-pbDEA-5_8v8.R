
fun <- \(x) {
    pb <- aggregateData(x)
    set.seed(1234)
    res <- bplapply(1:5, 
                      BPPARAM=MulticoreParam(3, RNGseed=1234, progress=FALSE), 
                      \(s){
        sw <- c(sample(which(pb$group_id=="A"),8),sample(which(pb$group_id=="B"),8))
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
                d
            }, error=\(e){ warning(e); NULL })
        })
        dplyr::bind_rows(d[!sapply(d, is.null)], .id="celltype")
    })
}
