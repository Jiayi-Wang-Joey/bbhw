suppressPackageStartupMessages({
    library(data.table)
})
fun <- \(res, truth, th) {
    
    getThStats <- function(ssres, truth, th=0.1, scores=NULL){
        m <- dplyr::bind_rows(ssres, .id="seed")
        if(is.null(scores))
            scores <- grep("^FDR|^padj\\.", colnames(m), value=TRUE)
        if(is.null(truth$celltype)) truth$celltype <- truth$cluster_id
        m <- merge(m, truth[,c("celltype","gene","isDEG")], 
                   by=c("celltype","gene"), all.x=TRUE)
        m <- m[which(!is.na(m$isDEG)),]
        m <- m[order(m$celltype, m$PValue),]
        scores <- intersect(scores, colnames(ssres[[1]]))
        ss <- dplyr::bind_rows(lapply(setNames(scores,scores), FUN=function(x){
            dplyr::bind_rows(lapply(split(m, m$celltype), \(m){
                m2 <- m[which(m[[x]]<th),]
                ret <- data.frame(TP=sum(m2$isDEG), FP=sum(!m2$isDEG))
                ret$precision=as.numeric(ret[1]/sum(ret))
                ret$recall=as.numeric(ret[1]/sum(m$isDEG))
                ret$calledSig <- (ret$TP + ret$FP)/length(ssres)
                ret$F1 <- 2/(1/ret$precision+1/ret$recall)
                ret
            }), .id="celltype")
        }), .id="method")
        ss <- reshape2::melt(ss, id.vars=c("method", "celltype", "TP", "FP",
                                           "calledSig","F1"),
                             measure.vars=c("precision", "recall", "F1"))
        ss2 <- ss[ss$variable=="F1",c("method","value")]
        ss2 <- aggregate(ss2[,2], ss2[,1,drop=FALSE], FUN=mean)
        ss2$x[is.na(ss2$x)] <- 0
        ss$method <- factor(ss$method, ss2$method[order(ss2$x)])
        ss
    }
    ss <- getThStats(res, truth, th)
    data.frame(ss, threshold=th)
}
