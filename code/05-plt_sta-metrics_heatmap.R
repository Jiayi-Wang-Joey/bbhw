suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(viridis)
    library(patchwork)
    library(ComplexHeatmap)
})

res <- lapply(args[[1]], readRDS)
df <- do.call(rbind,res)
setDT(df)
df <- df[threshold==0.05 & size=="12_2v2",]
df <- df[!grepl("raw|^FDR$", df$method),]
.plot <- \(df, type = "mean") {
    d <- unique(df[,.(method, celltype, variable, value)])
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
    d <- d[grep(type, d$variable),]
    d[,variable:= gsub(paste0(type,"."), "", variable, fixed=TRUE)]
    d[,variable := factor(variable, c("precision","recall", "F1"))]
    d$method <- renameScores(gsub("padj\\.", "", d$method))
    wPrior <- ifelse(grepl("^FDR\\.glb$|^FDR\\.loc$",d$method), "no prior", 
                     "using bulk prior")
    
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
    wPrior <- ifelse(grepl("^FDR\\.glb$|^FDR\\.loc$", mo), "no prior", "using bulk prior")
    cols <- setNames(list(viridis::magma(100), viridis::inferno(100), viridis::viridis(100)),
                     names(ml))
    
    hl <- lapply(names(ml), function(x) {
        m <- ml[[x]][mo, ]
        if (grepl("Precision", x, ignore.case = TRUE)) {
            col_fun <- circlize::colorRamp2(seq(0, 1, length.out=100), cols[[x]])
        } else {
            if (min(m, na.rm=TRUE) == max(m, na.rm=TRUE)) {
                col_fun <- circlize::colorRamp2(seq(0, 1, length.out=100), cols[[x]])
            } else {
                col_fun <- circlize::colorRamp2(
                    breaks = seq(min(m, na.rm=TRUE), max(m, na.rm=TRUE), length.out=100),
                    colors = cols[[x]]
                )
            }
        }
        rownames(m) <- ifelse(
            rownames(m) %in% c("FDR.glb", "FDR.loc"),
            paste0(rownames(m), " (no prior)"),
            rownames(m)
        )
        Heatmap(
            m, name = x,
            cluster_columns = FALSE, cluster_rows = FALSE,
            col = col_fun,
            column_title = x,
            column_names_gp = gpar(fontsize=9),
            row_split = wPrior,
            row_names_gp = gpar(fontsize=9)
        )
    })
    
    
    return(draw(Reduce("+", hl), merge=TRUE))
}
ps <- lapply(split(df, df$sim), \(fd) .plot(fd))


pdf(args[[2]], width = 7, height = 6)
for (p in ps) print(p)
dev.off()
