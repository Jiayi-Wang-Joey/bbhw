suppressPackageStartupMessages({
    library(muscat)
    library(BiocParallel)
    library(dplyr)
    library(SingleCellExperiment)
    library(scuttle)
})

pbDEAs <- readRDS(args$pb)
bulkDEA <- readRDS(args$bulk)
sce <- readRDS(args$sce)
pb <- aggregateData(sce)
rs <- sapply(assays(pb), rowSums)
props <- (1L+rs)/(rowSums(1L+rs))
rm(sce)
loc <- ifelse(wcs$loc=="loc", TRUE, FALSE)
bin <- wcs$bin
cor <- wcs$cor

lst <- lapply(pbDEAs, \(pbDEA) {
    pbDEA <- pbDEA[,c("celltype", "gene", "PValue", "logFC", "FDR")]
    pbDEA$FDR.global <- p.adjust(pbDEA$PValue, method="fdr")
    pbDEA <- dplyr::bind_rows(lapply(split(pbDEA, pbDEA$celltype), \(x){
        x$FDR.local <- p.adjust(x$PValue, method="fdr")
        x
    }))
    gi <- intersect(unique(pbDEA$gene), row.names(pb))
    pbDEA <- pbDEA[which(pbDEA$gene %in% gi),]
    pbDEA$ID <- paste0("H",seq_len(nrow(pbDEA)))
    name <- paste("padj", bin, gsub("gBH\\.","",cor),
                  wcs$loc,sep=".")
    res <- tryCatch(
        bbhw(pbDEA, bulkDEA, 
             pb=props, 
             bin.method=bin,
             local=loc,
             correction.method=cor),
        error=function(e){
            message(paste(name, " failed"))
            return(NULL)
    })

    if (!is.null(res)) {
        rownames(res) <- res$ID
        pbDEA[[name]] <- res[pbDEA$ID, "padj"]
    }
    if(bin == "sig") {
        res <- bbhw(pbDEA, bulkDEA, bin.method=bin,
                    local=loc, useSign = FALSE,
                    correction.method=cor)
        row.names(res) <- res$ID
        name <- paste("padj", bin, gsub("gBH\\.","",cor),
                       wcs$loc, "raw", sep=".")
        pbDEA[[name]] <- res[pbDEA$ID, "padj"]
    }

    data.frame(pbDEA, wcs)
})



saveRDS(lst, args$res)
