suppressPackageStartupMessages({
    library(splatter)
    library(scater)
    library(scran)
    library(dplyr)
})

set.seed(seed <- 20250717)


fun <- \() {
    
    ss <- c(1.2, 1.2, 0.5, 0.5, 0.8, 0)
    ts <- rep(1, 6)
    ci <- c("A1","A2","B1", "B2", "C", "D")
    bcs <- c(200,100,200,100,100,200)

    xs <- lapply(seq_len(6), \(i) {
        c <- ci[[i]]
        bc <- bcs[[i]]
        s <- ss[i]
        t <- ts[i]
        p <- newSplatPopParams(
            de.facLoc=t,
            de.prob=0.1,
            de.facScale=t,
            cde.prob=0.1,
            cde.facLoc=s,
            cde.facScale=s,
            group.prob=1,
            batchCells=bc,
            bcv.common=0.3,
            similarity.scale = 1,
            condition.prob=c(0.5, 0.5)
        )
        
        vcf <- mockVCF(n.samples=30, seed=seed)
        gff <- mockGFF(n.genes=6e3, seed=seed)
        
        x <- splatPopSimulate(
            params=p, vcf=vcf, gff=gff,
            verbose=FALSE, sparsify=FALSE)
        if(length(c)==1) {
            x$Group <- c
        } else {
            x$Group <- ifelse(x$Group == "Group1", c[1], c[2])
        }
        #x
        tab <- table(x$Sample, x$Condition)
        
        samples1 <- rownames(tab)[tab[, "Condition1"] != 0]
        samples2 <- rownames(tab)[tab[, "Condition2"] != 0]
        
        lookup1 <- setNames(paste0("sample_", sprintf("%02d", seq_along(samples1))), samples1)
        lookup2 <- setNames(paste0("sample_", sprintf("%02d", length(samples1) + seq_along(samples2))), samples2)
        
        x$sample_id <- x$Sample
        x$sample_id[x$Sample %in% samples1] <- lookup1[x$Sample[x$Sample %in% samples1]]
        x$sample_id[x$Sample %in% samples2] <- lookup2[x$Sample[x$Sample %in% samples2]]
        
        
        x
    })
    
    
    rds <- lapply(xs, rowData)
    # save ground truth
    names(rds) <- ci
    de <- lapply(seq_len(length(rds)), \(id) {
        x <- names(rds)[id]
        y <- rds[[id]]
        i <- grep("^ConditionDE", names(y))[1]
        j <- grep("^ConditionDE", names(y))[2]
        lfc <- abs(abs(log2(y[[j]]/y[[i]])))
        data.frame(gene=rownames(y),
                   celltype=x,
                   logFC=lfc)
    }) |> do.call(what=rbind)
    
    
    xs <- lapply(xs, \(.) {rowData(.) <- NULL; .})
    x <- do.call(cbind,xs)
    rm(xs)
    # standardize cell metadata
    colData(x) <- DataFrame(
        cluster_id=x$Group,
        sample_id=x$sample_id,
        Sample = x$Sample,
        group_id=x$Condition)
    for (. in names(colData(x)))
        x[[.]] <- factor(x[[.]])
    
    
    y <- counts(x)
    gs <- rowSums(y > 1) >= 10
    cs <- colSums(y > 0) >= 10
    x <- x[gs, cs]

    x <- logNormCounts(x)
    tbl <- modelGeneVar(x, block=x$sample_id)
    hvg <- (rowData(x)$bio <- tbl$bio ) > 0
    x <- runPCA(x, subset_row=hvg, ncomponents=10)
    x <- runUMAP(x, dimred= "PCA")
    x$group_id <- ifelse(x$group_id=="Condition1", "A", "B")
    return(list(sce=x, truth=de))

}



