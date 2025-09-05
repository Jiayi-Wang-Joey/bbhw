suppressPackageStartupMessages({
    library(muscat)
    library(SingleCellExperiment)
    library(BiocParallel)
    library(data.table)
    library(HDF5Array)
    library(tidyr)
})


set.seed(123)
fun <- \() {
    ref <- readRDS("data/ref/SCE.subset.annotated.rds")
    ref <- as(ref, "SingleCellExperiment")
    colData(ref) <- DataFrame(
        cluster_id=ref$celltype,
        sample_id=ref$sample)
    ref <- ref[,grepl("H", ref$sample_id)] 
    cd <- as.data.table(colData(ref))
    cd[,idx:=.I]
    sc <- unique(cd[,c("sample_id", "cluster_id")])
    idx <- sapply(seq_len(nrow(sc)), \(i) {
        dc <- cd[sample_id==sc[i,]$sample_id & 
                     cluster_id==sc[i,]$cluster_id]
        sample(dc$idx, min(300,nrow(dc)))
    })
    # subset sce
    ref <- ref[, unlist(idx)]
    ref <- prepSim(ref)
    
    probs <- list(cluster=c(0.8,0.2),
                  sample=rep(1/15,15),
                  group=c(0.5,0.5))
    print("sim1 and 2")
    sim <- bplapply(list(c(0.5,0.5), c(0.15,0.15)), #0.5, 0.1
                    BPPARAM=MulticoreParam(2, RNGseed=123), \(x){
                        simData(ref, nc=16000, ns=15, nk=2, 
                                probs=probs, p_type=0.2, force=TRUE,
                                dd=TRUE, p_dd=c(0.75,0,0.25,0,0,0), # prob 0.25
                                paired=TRUE, rel_lfc=x)
                    })
    sim1 <- sim[[1]] #A1, A2
    sim2 <- sim[[2]] #B1, B2
    print("sim3")
    sim3 <- simData(ref, nc=3400, ns=15, nk=1, 
                    probs=c(list(1),probs[-1]), force=TRUE,
                    dd=TRUE, p_dd=c(0.8,0,0.2,0,0,0), paired=TRUE, rel_lfc=0.3) # C lfc=0.3, prob=0.2
    print("sim4")
    sim4 <- simData(ref, nc=8000, ns=15, nk=1, 
                    probs=c(list(1),probs[-1]), force=TRUE,
                    dd=TRUE, p_dd=c(1,0,0,0,0,0), paired=TRUE, rel_lfc=0) # D
    rm(ref)
    gi1 <- metadata(sim1)$gene_info
    levels(gi1$cluster_id) <- levels(sim1$cluster_id) <- 
        c("A1", "A2")
    gi2 <- metadata(sim2)$gene_info
    levels(gi2$cluster_id) <- levels(sim2$cluster_id) <- 
        c("B1", "B2")
    gi3 <- metadata(sim3)$gene_info
    levels(gi3$cluster_id) <- levels(sim3$cluster_id) <- 
        c("C")
    levels(sim4$cluster_id) <- "D"
    rowData(sim1) <- rowData(sim2) <- rowData(sim3) <- rowData(sim4) <- NULL
    sim <- cbind(sim1,sim2,sim3,sim4)
    de <- rbind(gi1[which(gi1$category=="de"),], gi2[which(gi2$category=="de"),],
                gi3[which(gi3$category=="de"),])
    de <- de[,c("gene","cluster_id","logFC")]
    de <- de[which(de$logFC!=0),]
    colnames(de)[2] <- "celltype"
    
    # DR
    sim <- scater::logNormCounts(sim)
    sim <- scater::runPCA(sim, 10, ntop=3000)
    sim <- scater::runUMAP(sim, dimred="PCA", 
                           BNPARAM=BiocNeighbors::AnnoyParam())
    
    return(list(sce=sim, truth=de))
}
