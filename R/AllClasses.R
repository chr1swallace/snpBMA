setClass("snpBMAdata",
         representation(X="matrix",
                        Y="numeric",
                        family="character",
                        tags="character"),
         validity=function(object) {
           if(length(object@Y) != nrow(object@X))
             stop("phenotype vector Y must have length == nrow(X)")
         })
##' 
setClass("snpBMA",
         representation(nsnps="numeric",
                        snps="character",
                        groups="list",
                        bf="matrix",
                        models="Matrix"),
         validity=function(object) {
           if(nrow(object@bf) != nrow(object@models))
             stop("bf must have same nrow as models")
           if(length(unique(object@snps)) != ncol(object@models))
             stop("models must have ncol == length(unique(snps))")
         })
setMethod("show", signature="snpBMA",
          function(object) {
            nmod <- nrow(object@models)
            nsnp <- ncol(object@models)
            n.use <- apply(object@models,1,sum)
            cat("Bayes factors for",nmod,"models from",nsnp,"SNPs.\nSNPs per model:\n")
            print(table(n.use))
            cat("Top models:\n")
            print(top.models(object))
          })

top.models <- function(object,nmodels=6) {
  o <- head(order(object@bf[,2],decreasing=TRUE),nmodels)
  tm <- object@models[o,]
  col.drop <- which(apply(tm,2,sum)==0)
  if(length(col.drop))
    tm <- tm[,-col.drop]
  bf <- object@bf[o,]
  colnames(bf) <- paste("twologB10-phi",1:ncol(bf),sep="")
  cbind(as(tm,"matrix"),bf)
}

snp.summary.1 <- function(object) {
  wh.snp <- colnames(object@models)[apply(object@models==1,1,which)]
  bf <- object@bf
  rownames(bf) <- wh.snp
  colnames(bf) <- paste("twologB10-phi",1:ncol(bf),sep="")
  o <- order(bf[,2],decreasing=TRUE)
  bf[o,]
}

snp.subset <- function(object,snp,value) {
  j=which(colnames(object@models)==snp)
  i <- which(object@models[,j]==value)
  data.frame(model=apply(object@models[ i, -j],1,paste,collapse=""),
             bf=object@bf[i, 2])  
}

snp.summary.2 <- function(object.new,object.old) {
  for(snp in object.new@snps) {
    summ.old <- snp.subset(object.old,snp,0)
    summ.new <- snp.subset(object.new,snp,1)
    summ <- merge(summ.old,summ.new,by="model",suffixes=c(".old",".new"))
    mean(summ[,"bf.new"] - summ[,"bf.old"])
  }
  
}
