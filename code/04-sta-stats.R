fun <- \(ssres, truth) {

    getStats <- function(sl, truth, celltype=rep(1L, length(truth)), roundNd=NULL, noRankAt1=TRUE){
      dplyr::bind_rows(lapply(split(seq_along(truth), celltype), FUN=\(i){
        sl <- sl[i,,drop=FALSE]
        truth <- truth[i]
        if(sum(truth)<3) return(data.frame(nominal=numeric(0), recall=numeric(0), fdr=numeric(0)))
        dplyr::bind_rows(lapply(as.data.frame(sl), \(x){
          o <- order(x)
          if(noRankAt1) o <- intersect(order(x),which(x<1))
          d <- data.frame(nominal=c(x[o],1), label=c(truth[o],NA),
                          recall=c(cumsum(truth[o])/sum(truth),1),
                          fdr=c(cumsum(!truth[o])/seq_along(x[o]), 1-sum(truth)/length(truth)))
          if(!is.null(roundNd)){
            w <- which(d$nominal>0.25)
            d$recall <- round(d$recall, roundNd)
            d$fdr <- round(d$fdr, roundNd)
            d <- d[rev(seq_len(nrow(d))),]
            d <- d[!duplicated(d[,3:4]),]
          }
          d
        }), .id="score")
      }), .id="celltype")
    }
    m <- dplyr::bind_rows(ssres, .id="seed")
    m <- merge(m, truth, by=c("celltype","gene"), all.x=TRUE)
    m$isDEG <- !is.na(m$logFC.y)
    scores <- grep("^FDR|^padj\\.", colnames(m), value=TRUE)
    m <- m[order(m$celltype, m$PValue),]
    sts <- getStats(m[,intersect(scores,colnames(m))], m$isDEG, 
                    celltype=m$celltype, roundNd = NULL, noRankAt1=TRUE)
}


