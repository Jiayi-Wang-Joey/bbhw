args <- list(list.files("outs/sim/sta/", "^F1.*", full.names = TRUE), "plts/sim/das-p_hm.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(viridis)
})

res <- lapply(args[[1]], readRDS)
df <- do.call(rbind,res)
setDT(df)
df <- df[threshold==0.05  & sim=="muscat_LPS",]
df <- df[!grepl("raw|^FDR$", df$method),]
d <- unique(df[,.(method, celltype, TP, FP, 
                  calledSig, F1, variable, 
                  value, size, bin, loc, cor)])
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

d <- d[variable!="precision"]
d[,size:=factor(size)]
ggplot(d, aes(x = size, y = value, 
              group = method, col=method)) +
    geom_line() + 
    geom_point() +
    facet_grid2(variable ~ celltype, scales = "free", independent = "y") +
    scale_color_brewer(palette = "Paired") +
    theme_bw() +
    theme(panel.grid = element_blank())


