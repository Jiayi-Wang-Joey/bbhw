suppressPackageStartupMessages({
    library(muscat)
    library(SingleCellExperiment)
    library(BiocParallel)
    library(HDF5Array)
})

# ss <- c(1.2, 1.2, 0.5, 0.5, 0.8, 0)
# ts <- rep(1, 6)
# ci <- c("A1","A2","B1", "B2", "C", "D")
# bcs <- c(200,100,200,100,100,200)
set.seed(123)
fun <- \() {
    ref <- readRDS("data/ref/PBMC.subset.annotated.se.rds")
    assay(ref, "counts", withDimnames = FALSE) <- HDF5Array(
        filepath = "data/ref/SCE.subset.annotated..assays.h5",
        name = "assay001"
    )
    ref <- as(ref, "SingleCellExperiment")
    colData(ref) <- DataFrame(
        cluster_id=ref$celltype,
        sample_id=ref$sample)
    
    ref <- prepSim(ref)
    probs <- list(cluster=c(0.8,0.2),
                  sample=rep(0.25,15),
                  group=c(0.5,0.5))
    
    sim <- bplapply(list(c(0.3,0.3), c(0.05,0.05)), 
                    BPPARAM=MulticoreParam(2, RNGseed=123), \(x){
                        simData(ref, nc=12000, ns=15, nk=2, 
                                probs=probs, p_type=0.2, force=TRUE,
                                dd=TRUE, p_dd=c(0.80,0,0.2,0,0,0), 
                                paired=TRUE, rel_lfc=x)
                    })
    sim1 <- sim[[1]] #A1, A2
    sim2 <- sim[[2]] #B1, B2
    sim3 <- simData(ref, nc=2500, ns=15, nk=1, 
                    probs=c(list(1),probs[-1]), force=TRUE,
                    dd=TRUE, p_dd=c(0.9,0,0.1,0,0,0), paired=TRUE, rel_lfc=0.15) # C
    sim4 <- simData(ref, nc=1200, ns=15, nk=1, 
                    probs=c(list(1),probs[-1]), force=TRUE,
                    dd=TRUE, p_dd=c(1,0,0,0,0,0), paired=TRUE, rel_lfc=0) # D
    
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