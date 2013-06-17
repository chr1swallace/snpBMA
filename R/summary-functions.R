##' Show top models in a snpBMA object
##'
##' @title snpBMA
##' @param object Object of class snpBMA
##' @param nmodels number of models to show
##' @return matrix of top models and their Bayes Factors
##' @author Chris Wallace
top.snpBMA <- function(object,nmodels=6L, keep.snps=FALSE) {
  if(is.na(nmodels) || nmodels>nrow(object@bf))
    nmodels <- nrow(object@bf)
  o <- head(order(object@bf[,2],decreasing=TRUE),nmodels)
  tm <- object@models[o,]
  if(!keep.snps) {
    col.drop <- which(apply(tm,2,sum)==0)
    if(length(col.drop))
      tm <- tm[,-col.drop]
  }
  bf <- object@bf[o,]
  colnames(bf) <- paste("twologB10-phi",1:ncol(bf),sep="")
  cbind(as(tm,"matrix"),bf)
}

post.snpBMA <- function(object, prior) {
  onames <- model.names(object@models)
  obf <- object@bf
  snps.all <- snps0(object)
  snps <- lapply(strsplit(onames, "-"), function(x) snps.all[x])
  rownames(obf) <- lapply(snps, paste, collapse="-")
  pp <- exp(obf/2) * prior
  colnames(pp) <- paste("PP_phi=phi.", 1:ncol(pp), sep="")
  return(pp)
}

##' Show top models in a snpBMA object
##'
##' @title snpBMA
##' @param object Object of class snpBMA
##' @param nmodels number of models to show
##' @return matrix of top models and their Bayes Factors
##' @author Chris Wallace
top.snpBMAlist <- function(object,priors,nmodels=6L) {
  n <- length(object)
  priors2 <- priors[1:n]/sapply(as.list(1:n), function(i) { object[[i]]@nmodels })
  pp <- lapply(1:n, function(i) {
    post.snpBMA(object[[i]], priors2[i])
  })
  cn <- colnames(pp[[1]])
  pp <- do.call("rbind",pp)
  pp <- t(t(pp) / colSums(pp))
  pp <- pp[order(pp[,2], decreasing=TRUE),]
  colnames(pp) <- cn
  head(x=pp, n=nmodels)
}

top.snps <- function(object, ...) {
  ts <- top.models(object, ...)
  unique(unlist(strsplit(rownames(ts), "-")))
}

##' Summary Bayes Factors for each SNP
##'
##' The summary Bayes Factor is the mean of the Bayes Factor for all
##' models containing that SNP, which is valid if each model has equal
##' prior
##' @title Summary Bayes Factors for each SNP
##' @param object Object of class snpBMA
##' @param do.order if TRUE (the default), reorder the returned summary of snp Bayes Factors, so that those with the strongest support are shown first
##' @return a Matrix of class \code{\linkS4class{dgeMatrix}} with each row corresponding to a SNP, and each column to a different value of phi
##' @author Chris Wallace
##' @export
snp.summary <- function(object, do.order=TRUE) {
  snp.bf <- Matrix::t(object@models) %*% object@bf / Matrix::colSums(object@models)
  if(do.order) 
    snp.bf <- snp.bf[ order(snp.bf[,2], decreasing=TRUE), ]
  colnames(snp.bf) <- paste("twologB10-phi",1:ncol(snp.bf),sep="")
  return(snp.bf)
}

snp.subset <- function(object,model.str, snp,value) {
  j=which(colnames(object@models)==snp)
  nc <- ncol(object@models)
  i <- which(object@models[,j]==value)
  m.left <- m.right <- NULL
  if(j>1)
    m.left <- substr(model.str[i],1,j-1)
  if(j<nc)
    m.right <- substr(model.str[i],j+1,nc)
  ret <- object@bf[i, ]
  rownames(ret) <- paste(m.left,m.right,sep="")
  return(ret)
  ## data.frame(model=paste(m.left,m.right,sep=""), #apply(object@models[ i, -j],1,paste,collapse=""),
  ##            bf=)  
}

snp.summary.2 <- function(object.new,object.old, do.order=TRUE) {
  ## check there is a one SNP difference in models
  snp.diff <- unique(object.new@nsnps) - unique(object.old@nsnps)
  if(snp.diff != 1)
    stop("object.new should have 1 more SNP in all models than object.old.")
  snps <- unique(object.new@snps)
  bf.snp <- matrix(NA,length(snps), ncol(object.new@bf))
  rownames(bf.snp) <- snps
  nc <- ncol(bf.snp)
  models.old <- apply(object.old@models,1,paste,collapse="")
  models.new <- apply(object.new@models,1,paste,collapse="")
  for(snp in object.new@snps) {
    cat(".")
    summ.old <- snp.subset(object.old,models.old,snp,0)
    summ.new <- snp.subset(object.new,models.new,snp,1)
    rn <- intersect(rownames(summ.old),rownames(summ.new))
    summ <- cbind(summ.old[rn,], summ.new[rn,])
    bf.snp[snp, ] <-  colMeans(summ[,-c(1:nc)] - summ[,c(1:nc)])
  }  
  if(do.order) 
    bf.snp <- bf.snp[ order(bf.snp[,2], decreasing=TRUE), ]
}
