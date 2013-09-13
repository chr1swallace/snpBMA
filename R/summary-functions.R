##' Show top models in a snpBMA object
##'
##' @title snpBMA
##' @param object Object of class snpBMA
##' @inheritParams top.snpBMAlist
##' @return matrix of top models and their Bayes Factors
##' @author Chris Wallace
top.snpBMA <- function(object,nmodels=6L, keep.snps=FALSE, what="models") {
##   if(what=="snps") {
##     ## convert to bf for snps, assume priors all equal
##     bf <- snp.summary(object, do.order=FALSE)
##     tm <- as(rownames(bf),"matrix")
##   } else {
    bf <- object@bf
    tm <- object@models
##   }
  if(is.na(nmodels) || nmodels>nrow(bf))
    nmodels <- nrow(bf)
  o <- head(order(bf[,2,drop=FALSE],decreasing=TRUE),nmodels)
  tm <- tm[o,,drop=FALSE]
  if(what=="models" & !keep.snps) {
    col.drop <- which(colSums(tm)==0)
    if(length(col.drop))
      tm <- tm[,-col.drop,drop=FALSE]
  }
  bf <- bf[o,,drop=FALSE]
  colnames(bf) <- paste("twologB10-phi",1:ncol(bf),sep="")
  if(what=="snps") {
    rownames(bf) <- tm[,1]
    return(bf)
  }
  return(cbind(as(tm,"matrix"),as(bf,"matrix")))
}
##' post.snpBMA
##'
##' Calculate posterior probs given BF and prior odds
##' @title Internal, post.snpBMA
##' @param object Object of class snpBMA
##' @param priors prior odds for models containing a fixed number of SNPs
##' @param snpsep character string used to separate SNP names when naming models
##' @return pp matrix
##' @author Chris Wallace
post.snpBMA <- function(object, priors, snpsep="/") {
  onames <- model.names(object@models)
  obf <- object@bf
  snps.all <- snps0(object)
  snps <- lapply(strsplit(onames, "-"), function(x) snps.all[x])
  pp <- exp(obf/2) * priors
  colnames(pp) <- paste("PP_phi=phi.", 1:ncol(pp), sep="")
  return(new("modelSummary",.Data=pp, snps=snps, snpsep=snpsep))
}

## tag.expand <- function(tags, tags.99) {

  
## }

##' Internal function, logsum
##'
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##' @title logsum (copied from coloc package)
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
##' @author Claudia Giambartolomei
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}
##' Internal function, logdiff
##'
##' This function calculates the log of the difference of the exponentiated
##' logs taking out the max, i.e. insuring that the difference is not negative
##' @title logdiff (copied from coloc package)
##' @param x numeric
##' @param y numeric
##' @return max(x) + log(exp(x - max(x,y)) - exp(y-max(x,y)))
##' @author Chris Wallace
logdiff <- function(x,y) {
  my.max <- pmax(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}


##' Show top models in a snpBMAlist object
##'
##' @title top.snpBMAlist
##' @param object Object of class snpBMAlist
##' @param priors vector of priors from which prior odds for models containing a fixed number of SNPs may be derived
##' @param nmodels number of models to show
##' @param min.pp ("minimim posterior probability") optional, if present, limit to models with psterior probabilities > min.pp
##' @param min.cpp ("minimum cumulative posterior probability") optional, if present, limit to the best models which together account for at least min.cpp of the posterior probability
##' @param what \code{models} or \code{snps}, indicating whether to output the top models or SNPs
##' @param ... additional arguments passed to \link{post.snpBMA}
##' @return matrix of top models and their Bayes Factors
##' @author Chris Wallace
top.snpBMAlist <- function(object,priors,nmodels=6L, min.pp=NULL, min.cpp=NULL, what="models", tags=NULL, groups=NULL, ...) {
  n <- length(object)

  if(!is.null(tags)) {
#    tags <- tags[ tags %in% object[[1]]@snps ]
    groups <- split(names(tags),tags)
  }

  ##   priors2 <- priors[1:n]
  ## if(is.null(prior.odds) & is.null(prior.prob))
  ##   stop("Must specify prior odds OR prior probs")
  ## if(!is.null(prior.prob)) {
  ##   priors2 <- priors[1:n]/sapply(as.list(1:n), function(i) { bma.list[[i]]@nmodels })
  ## } else {
  ##   priors2 <- prior.odds
  ## }
  pp <- lapply(1:n, function(i) {
    post.snpBMA(object=object[[i]], priors=priors[i], ...)
  })
  if(what=="snps") {
    r <- unlist(lapply(pp,nrow))
  }
  pp <- do.call("stack",pp)
  snps <- pp@snps
  snpsep <- pp@snpsep
  pp <- t(t(pp) / colSums(pp))
  if(what=="snps") {
    snps <- unique(unlist(lapply(object, function(o) o@snps)))
    models <- Matrix(0,sum(r),length(snps),dimnames=list(NULL,snps))
    rc <- c(0,cumsum(r))
    for(i in 1:n) {
      r <- (rc[i]+1):rc[i+1]
      m <- object[[i]]@models
      models[r,colnames(m)] <- m
    }
    if(!is.null(groups)) {
      newmodels <- Matrix(0,nrow(models),length(groups),dimnames=list(NULL,names(groups)))
      for(g in names(groups))
        newmodels[,g] <- rowSums(models[,groups[[g]],drop=FALSE])
      models <- newmodels
    }
    pp <- Matrix::t(models) %*% pp
    snps <- colnames(models)    
    pp <- new("modelSummary", .Data=pp, snps=as.list(snps), snpsep=snpsep)
 }
  o <- order(pp[,2],decreasing=TRUE)
  pp <- pp[o,,drop=FALSE]
 
  ## return values
  if(!is.null(min.pp))
    return(pp[ pp[,2]>=min.pp, ])
  if(!is.null(min.cpp)) {
    cpp <- cumsum(pp[,2])
    wh <- which(cpp>=min.cpp)[1]
    return(pp[ 1:wh, ])
  }      
  if(is.na(nmodels) || nmodels==0)
    return(pp)
  return(pp[1:nmodels,])
}


##' Extract SNPs from top.models
##'
##' @title top.snps
##' @param object snpBMA or snpBMAlist object
##' @param ... additional arguments to pass to top.models, notably \code{nmodels}: the number of models to include
##' @return character vector of SNP names
##' @author Chris Wallace
##' @export
top.snps <- function(object, snpsep="/", ...) {
  ts <- top.models(object, snpsep=snpsep, ...)
  unique(unlist(strsplit(rownames(ts), snpsep)))
}

##' Summary Bayes Factors for each SNP
##'
##' The summary Bayes Factor is the sum of the Bayes Factor for all
##' models containing that SNP, which is valid if each model has equal
##' prior
##' @title Summary Bayes Factors for each SNP
##' @param object Object of class snpBMA
##' @param do.order if TRUE (the default), reorder the returned summary of snp Bayes Factors, so that those with the strongest support are shown first
##' @return a Matrix of class \code{\linkS4class{dgeMatrix}} with each row corresponding to a SNP, and each column to a different value of phi
##' @author Chris Wallace
##' @export
snp.summary <- function(object, do.order=TRUE) {
  if(is(object,"snpBMAstrat")) {
    bf.list <- lapply(object@.Data, snp.summary, do.order=FALSE)
    snp.bf <- bf.list[[1]]
    for(i in 2:length(bf.list))
      snp.bf <- snp.bf + bf.list[[i]]
    if(do.order)
      snp.bf <- snp.bf[ order(snp.bf[,2], decreasing=TRUE), ]
    return(snp.bf)
  }
  snp.bf <- Matrix::t(object@models) %*% object@bf
  ## add in tags
  toadd <- setdiff(names(object@tags),object@tags)
  if(length(toadd)) {
    bf.toadd <- snp.bf[ object@tags[toadd], ]
    rownames(bf.toadd) <- toadd
    snp.bf <- rbind(snp.bf, bf.toadd)
  }
  ## drop snps not in any models
  rs <- rowSums(abs(snp.bf))
  snp.bf <- snp.bf[rs>0,,drop=FALSE]
  if(do.order) 
    snp.bf <- snp.bf[ order(snp.bf[,2], decreasing=TRUE), ]
  colnames(snp.bf) <- paste("twologB10-phi",1:ncol(snp.bf),sep="")
  return(snp.bf)
}

nsnps.posterior <- function(object) {
  nsnps <- sapply(object@snps,length)
  maxsnps <- max(nsnps)
  t(sapply(1:maxsnps, function(i) {
    wh <- which(nsnps==i)
    colSums(object[wh,,drop=FALSE])
  })) ## posterior distribution on number of SNPs  
}

nsnps.expected <- function(object) {
  nsnps <- sapply(object@snps,length)
  t(nsnps) %*% object ## expected number of snps
}

snps.pp <- function(object, n=NULL, rescale=FALSE, byn=FALSE) {
  nsnps <- sapply(object@snps,length)
  if(byn && is.null(n)) { ## split by nsnps, do all of them
    maxsnps <- max(nsnps)
    PPisnp <- lapply(1:maxsnps, function(i) cbind(nsnps=i,as(snps.pp(object,i,rescale=rescale),"matrix")))
    return(do.call("rbind",PPisnp))
  }
  if(!is.null(n)) { ## split, do subset of snps
    wh <- which(nsnps==n)
    object <- object[wh,]
    nsnps <- nsnps[wh]
  }
  snps <- unique(unlist(object@snps))
  models <- Matrix(0,nrow(object),length(snps),dimnames=list(1:nrow(object),snps))
  index <- cbind(rep(1:length(nsnps), times=nsnps),
                 match(unlist(object@snps),colnames(models)))
  models[index] <- 1
  pp <- object@.Data
  if(rescale)
    pp <- t(t(pp)/colSums(pp))
  return(t(models) %*% pp)
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
