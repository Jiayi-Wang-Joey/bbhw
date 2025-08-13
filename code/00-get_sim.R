source(args[[1]])
res <- fun()

saveRDS(res$sce, args$sce)
saveRDS(res$truth, args$truth)