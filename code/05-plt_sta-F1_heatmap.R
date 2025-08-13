args <- list(list.files("outs/sim/sta/", "^F1.*", full.names = TRUE), "plts/sim/das-p_hm.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(viridis)
})

res <- lapply(args[[1]], readRDS)
df <- do.call(rbind,res)
setDT(df)
df <- df[threshold==0.05 & size=="5_2v2" & sim=="splatter",]
df <- df[!grepl("raw|^FDR$", df$method),]
d <- unique(df[,.(method, celltype, TP, FP, calledSig, F1, variable, value)])
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
d$wPrior <- ifelse(grepl("^FDR\\.glb$|^FDR\\.loc$",d$method), "no prior", 
                   "using bulk prior")
td <- d[celltype!="D"]
ggplot(td, aes(x = celltype, y = method, fill = value)) +
    geom_tile(color = "white") +
    facet_grid(wPrior ~ variable, scales="free") +
    scale_fill_gradientn(colors = viridis(100)) +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_rect(fill = "grey90", color = NA)
    )
    