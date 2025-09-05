suppressPackageStartupMessages({
    library(ggplot2)
    library(SingleCellExperiment)
    library(patchwork)
    library(ggrastr)
    library(RColorBrewer)
})
fs <- args[[1]][grepl("PBMC", args[[1]])]
sce <- readRDS(fs)
df <- data.frame(reducedDim(sce, "UMAP"), colData(sce))

aes <- list(
    geom_point_rast(aes(UMAP1, UMAP2), 
                    shape=16, alpha=0.1, size=0.5),
    guides(col=guide_legend(
        ncol=4, title.position="top",
        override.aes=list(alpha=1, size=2))),
    theme_minimal(6), theme(
        aspect.ratio=1,
        plot.margin=margin(),
        axis.text=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(hjust=0),
        panel.border=element_rect(fill=NA),
        legend.key.size=unit(0.25, "lines")))


nk <- length(unique(df$cluster_id))
gg <-  ggplot(df, aes(col=group_id)) +
    scale_color_manual(values=c("royalblue", "tomato")) +
    ggplot(df, aes(col=cluster_id)) +
    scale_color_manual(values=colorRampPalette(brewer.pal(12, "Dark2"))(nk)) +
    plot_layout(nrow=1, guides="collect") & 
    aes & theme(
        plot.margin=margin(0),
        #legend.position="bottom",
        legend.justification=c(0.5, 1),
        legend.box.spacing=unit(0, "pt"),
        #plot.tag=element_text(size=9, face="bold")
        ) 

ggsave(args[[2]], gg, height = 3, width =5)