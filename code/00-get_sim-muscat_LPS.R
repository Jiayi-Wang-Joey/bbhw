suppressPackageStartupMessages({
    library(muscat)
    library(SingleCellExperiment)
    library(BiocParallel)
})
set.seed(123)
fun <- \() {
    ref <- readRDS("data/ref/week13.SCE.processed.rds")
    ref <- ref[,ref$group_id=="WT"]
    ref$cluster_id <- ref$cluster2
    ref$group_id <- NULL
    ref <- prepSim(ref)
    probs <- list(cluster=c(0.8,0.2),
                  sample=rep(0.25,8),
                  group=c(0.5,0.5))
    
    sim <- bplapply(list(c(0.2,0), c(0,0.2)), 
                    BPPARAM=MulticoreParam(2, RNGseed=123), \(x){
        simData(ref, nc=16000, ns=8, nk=2, probs=probs, p_type=0.2, force=TRUE,
                dd=TRUE, p_dd=c(0.80,0,0.2,0,0,0), paired=TRUE, rel_lfc=x)
    })
    sim1 <- sim[[1]]
    sim2 <- sim[[2]]
    sim3 <- simData(ref, nc=3200, ns=8, nk=1, 
                    probs=c(list(1),probs[-1]), force=TRUE,
                    dd=TRUE, p_dd=c(0.9,0,0.1,0,0,0), paired=TRUE, rel_lfc=0.15)
    sim4 <- simData(ref, nc=12800, ns=8, nk=1, 
                    probs=c(list(1),probs[-1]), force=TRUE,
                    dd=TRUE, p_dd=c(1,0,0,0,0,0), paired=TRUE, rel_lfc=0)
    
    gi1 <- metadata(sim1)$gene_info
    levels(gi1$cluster_id) <- levels(sim1$cluster_id) <- 
        c("largeAffected", "smallUnaffected")
    gi2 <- metadata(sim2)$gene_info
    levels(gi2$cluster_id) <- levels(sim2$cluster_id) <- 
        c("largeUnaffected1", "smallAffected1")
    gi3 <- metadata(sim3)$gene_info
    levels(gi3$cluster_id) <- levels(sim3$cluster_id) <- 
        c("smallAffected2")
    levels(sim4$cluster_id) <- "largeUnaffected2"
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