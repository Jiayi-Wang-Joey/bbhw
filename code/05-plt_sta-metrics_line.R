#args <- list(list.files("outs/sim/sta/", "^metrics*", full.names = TRUE), "plts/sim/das-p_hm.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(viridis)
    library(ggh4x)
    library(RColorBrewer)
    library(patchwork)
})

#fs <- args[[1]][grepl("^metrics.*\\.rds$", basename(args[[1]]))]
res <- lapply(args[[1]], readRDS)
df <- do.call(rbind,res)
setDT(df)
df <- df[threshold==0.05  & size != "12_2v2"]
df <- df[!grepl("raw|^FDR$", df$method),]

.plot <- \(df, type = "mean") {
    d <- unique(df[,.(method, celltype, variable, 
                      value, size, bin, loc, cor)])
    d <- d[grep(type, d$variable),]
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
    d[,variable:= gsub(paste0(type,"."), "", variable, fixed=TRUE)]
    d[,variable := factor(variable, c("precision","recall", "F1"))]
    d <- d[variable!="precision"]
    d[,size:=factor(size)]
    
    aes <- list(
        facet_grid2(variable ~ celltype, scales = "free", independent = "y"),
        guides(col=guide_legend(
            ncol=4, title.position="top",
        )),
        theme_minimal(), theme(
            plot.margin=margin(),
            panel.grid=element_blank(),
            panel.border=element_rect(fill=NA),
            legend.key.size=unit(0.25, "lines"),
            axis.text.x = element_text(angle = 45, hjust = 1) 
        )
    )
    
    d <- d[method %in% c("FDR.loc", "PAS.LSL.loc", "sig.LSL.loc", "asNA.LSL.loc")]
    d[method=="FDR.loc", method:="FDR.loc (no prior)"]
    d[,method:=factor(method, levels = c("FDR.loc (no prior)", "PAS.LSL.loc", 
                                         "sig.LSL.loc", "asNA.LSL.loc"))]
    n <- length(unique(d$method))
    ggplot(d, aes(x = size, y = value, 
                  group = method, col=method)) +
        geom_line(alpha=0.8) + geom_point(size=1.5, alpha=0.8) +
        scale_color_manual(values=colorRampPalette(brewer.pal(12, "Paired"))(n)) +
        plot_layout(nrow=1, guides="collect") & 
        aes & theme(legend.position="bottom")
    
}


ps <- lapply(split(df, df$sim), \(fd) .plot(fd))

pdf(args[[2]], width = 10, height = 4)
for (p in ps) print(p)
dev.off()
