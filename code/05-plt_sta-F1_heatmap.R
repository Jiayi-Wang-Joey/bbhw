args <- list(list.files("outs/sim/sta/", "^F1.*", full.names = TRUE), "plts/sim/das-p_hm.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(viridis)
    library(ComplexHeatmap)
})

res <- lapply(args[[1]], readRDS)
df <- do.call(rbind,res)
setDT(df)
df <- df[threshold==0.05 & size=="12_2v2" & sim=="splatter",]
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
wPrior <- ifelse(grepl("^FDR\\.glb$|^FDR\\.loc$",d$method), "no prior", 
                 "using bulk prior")
d <- d[celltype!="D"]

ml <- lapply(split(d, d$variable), FUN=function(x){
    m <- reshape2::dcast(x, formula = celltype~method, value.var="value")
    row.names(m) <- m[,1]
    t(m[,-1])
})

m <- do.call(cbind, lapply(ml, \(x){
    rm <- matrixStats::colMedians(x, na.rm=TRUE)
    for(i in seq_along(rm)){
        x[which(is.na(x[,i])),i] <- rm[i]
    }
    x
}))

mo <- row.names(m)[order(rowMeans(m, na.rm=TRUE) + 
                             matrixStats::rowMedians(m, na.rm=TRUE))]
wPrior <- ifelse(grepl("^FDR\\.glb$|^FDR\\.loc$",mo), "no prior", "using bulk prior")
cols <- setNames(list(viridis::magma(100), viridis::inferno(100), viridis::viridis(100)),
                 names(ml))
hl <- lapply(names(ml), FUN=function(x){
    m <- ml[[x]][mo,]
    if (min(m, na.rm=TRUE) == max(m, na.rm=TRUE)) {
        col_fun <- function(z) rep(cols[[x]][100], length(z))
    } else {
        col_fun <- circlize::colorRamp2(
            breaks = seq(min(m, na.rm=TRUE), max(m, na.rm=TRUE), length.out=100),
            colors = cols[[x]]
        )
    }
    
    Heatmap(
        m, name=x,
        cluster_columns = FALSE, cluster_rows = FALSE,
        col = col_fun,               
        column_title = x,
        column_names_gp = gpar(fontsize=9),
        row_split = wPrior,
        row_names_gp = gpar(fontsize=9)
    )
})

draw(Reduce("+", hl), merge=TRUE)
