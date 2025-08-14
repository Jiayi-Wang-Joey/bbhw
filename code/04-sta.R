suppressPackageStartupMessages({
    library(data.table)
})

truth <- readRDS(args$truth)
res <- readRDS(args$res)
source(args$fun)
res$cluster_id <- res$celltype
if(wcs$sim=="splatter") truth$isDEG <- truth$logFC != 0
if(wcs$sim=="muscat_LPS") truth$isDEG <- truth$logFC != 0 & !is.na(truth$logFC)
th <- c(0.05,0.1,0.25)
sta <- lapply(th, \(x) fun(res, truth, x)) |> do.call(what=rbind)
sta <- data.frame(sta, wcs)
sta$score <- res$score[1]
saveRDS(sta, args$sta)