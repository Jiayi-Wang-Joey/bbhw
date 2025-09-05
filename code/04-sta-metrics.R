suppressPackageStartupMessages({
    library(data.table)
})
fun <-  function(ssres, truth=NULL, th=0.1, scores=NULL) {
    .fun <- \(m, truth, th, scores=NULL) {
        if(is.null(scores))
            scores <- grep("^FDR|^padj\\.", colnames(m), value=TRUE)
        if(is.null(m$isDEG)){
            m <- merge(m, truth[,c("celltype","gene","isDEG")], 
                       by=c("celltype","gene"))
        }
        
        m <- m[which(!is.na(m$isDEG)),]
        m <- m[order(m$celltype, m$PValue),]
        ss <- dplyr::bind_rows(lapply(setNames(scores,scores), FUN=function(x){
            dplyr::bind_rows(lapply(split(m, m$celltype), \(m){
                m2 <- m[which(m[[x]]<th),]
                ret <- data.frame(TP=sum(m2$isDEG), FP=sum(!m2$isDEG))
                ret$precision=as.numeric(ret[1]/sum(ret)) 
                ret$recall=as.numeric(ret[1]/sum(m$isDEG))
                ret$F1 <- 2/(1/ret$precision+1/ret$recall)
                ret
            }), .id="celltype")
        }), .id="method")
    }
    ss <- lapply(ssres, \(x) .fun(x, truth, th, scores=scores))
    dt <- rbindlist(ss, idcol = "seed")
    dt[is.na(F1), F1 := 0]
    td <- dt[, .(
        mean.precision = mean(precision, na.rm=TRUE),
        sem.precision = sd(precision, na.rm=TRUE)/sqrt(length(precision)),
        median.precision = median(precision, na.rm=TRUE),
        mad.precision = median(abs(precision-median(precision, na.rm=TRUE))),
        mean.recall = mean(recall),
        sem.recall = sd(recall)/sqrt(length(recall)),
        median.recall = median(recall),
        mad.recall = median(abs(recall-median(recall, na.rm=TRUE))),
        mean.F1 = mean(F1),
        sem.F1 = sd(F1)/sqrt(length(F1)),
        median.F1 = median(F1),
        mad.F1 = median(abs(F1-median(F1, na.rm=TRUE))),
        TP    = as.numeric(median(TP)),
        mean.FP    = as.numeric(mean(FP))
    ), by = .(method, celltype)]
    
    td1 <- melt(
        td,
        measure.vars = paste(rep(c("mean","median"), each=3), c("precision", "recall", "F1"), sep="."),  
        id.vars = c("method", "celltype")
    )
    td2 <- melt(
        td,
        measure.vars = paste(rep(c("sem","mad"), each=3), c("precision", "recall", "F1"), sep="."),  
        id.vars = c("method", "celltype"),
        value.name = "var"
    )
    td1$var <- td2$var
    td1 <- td1[!grepl("Unaffected", td1$celltype),]
    td1 <- td1[!grepl("D", td1$celltype),]
    data.table(td1, threshold=th)
}
