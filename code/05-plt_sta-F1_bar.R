args <- list(list.files("outs/sim/sta/", "^F1.*", full.names = TRUE), "plts/sim/das-p_hm.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(viridis)
})

res <- lapply(args[[1]], readRDS)
df <- do.call(rbind,res)
setDT(df)
df <- df[threshold==0.05 & size!="12_2v2" & sim=="splatter",]
df <- df[!grepl("raw|^FDR$", df$method),]
d <- unique(df[,.(method, celltype, TP, FP, calledSig, F1, variable, value, size)])
renameScores <- function(st, rmLoc=FALSE, rmRaw=TRUE, rmSig=FALSE){
    if(is.data.frame(st)){
        if(!is.null(st$score))
            st$score <- renameScores(st$score, rmLoc=rmLoc, rmRaw=rmRaw, rmSig=rmSig)
        if(!is.null(st$method))
            st$method <- renameScores(st$method, rmLoc=rmLoc, rmRaw=rmRaw, rmSig=rmSig)
        return(st)
    }
    score <- factor(st)
    levels(score) <- gsub("^padj\\.","",levels(score))
    levels(score) <- gsub("Global", "glb", levels(score), ignore.case = TRUE)
    levels(score) <- gsub("Local", "loc", levels(score), ignore.case = TRUE)
    if(rmLoc) levels(score) <- gsub("\\.loc|\\.glb","",levels(score))
    if(rmRaw) levels(score) <- gsub("\\.raw","",levels(score))
    if(rmSig) levels(score) <- gsub("sig\\.","",levels(score))
    score
}

d$method <- renameScores(gsub("padj\\.", "", d$method))
wPrior <- ifelse(grepl("^FDR\\.glb$|^FDR\\.loc$",d$method), "no prior", 
                 "using bulk prior")
d <- d[celltype!="D" & method %in% c("FDR.loc", "PAS.LSL.loc", "sig.IHW.glb")]
d[,size:=factor(size)]
ggplot(d, aes(method, value, fill = size)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    facet_grid(variable~celltype)
