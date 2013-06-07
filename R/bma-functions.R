library(BMA)
library(ggplot2)
library(pheatmap)
library(reshape)

##' Calculate D' and R^2 for a SnpMatrix object and display as a heatmap
##'
##' Side effect: displays plot on current graphics device
##' @title show.ld
##' @inheritParams make.groups
##' @return invisibly returns matrix with D' in upper.tri() entries and R^2 in lower.tri() entries
##' @author Chris Wallace
show.ld <- function(snps=colnames(X), samples=rownames(X), X=snp.data,
                    lines.limit=20) {
 my.colours <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
                 "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  if(!is(X,"SnpMatrix"))
    stop("X must be a SnpMatrix")
 breaks <- NULL
 groups <- NULL
  if(is.list(snps)) {
    breaks <- sapply(snps,length)
    groups <- rep(1:length(snps), times=breaks)
    snps <- unlist(snps)
  }
  LD <- ld(X[samples,snps],
           depth=length(snps)-1,
           stat=c("D.prime","R.squared"),
           symmetric=TRUE)
  ld <- as.matrix(LD$R.squared)
  ld[upper.tri(ld)] <- t(as.matrix(LD$D.prime))[ upper.tri(ld) ]
  diag(ld) <- 1
 if(!is.null(groups))
   colnames(ld) <- rownames(ld) <- paste(groups,colnames(ld),sep="/")
 
  df <- melt(as.matrix(ld))
  df$X1 <- factor(df$X1, levels=colnames(ld))
  df$X2 <- factor(df$X2, levels=colnames(ld))
 n <- max(as.numeric(df$X1))
 p <- ggplot(df) + geom_tile(mapping=aes(x=X1,y=X2,fill=value), linetype=1) +
    scale_fill_gradientn(colours = my.colours) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text = element_text(colour="black"),
            axis.title = element_blank())
 if(n<lines.limit) {
   p <- p + geom_hline(yintercept=(0:n)+0.5, col="grey60") +
     geom_vline(xintercept=(0:n)+0.5, col="grey60")
 }
 
if(!is.null(breaks)) {
  x <- c(0,cumsum(breaks),max(as.numeric(df$X1))) + 0.5
  n <- length(x)
  xf <- rbind(data.frame(x=x[-n], y=x[-n], xend=x[-1], yend=x[-n]),
              data.frame(x=x[-n], y=x[-n], xend=x[-n], yend=x[-1]),
              data.frame(x=x[-1], y=x[-1], xend=x[-1], yend=x[-n]),
              data.frame(x=x[-1], y=x[-1], xend=x[-n], yend=x[-1]))
              
  show(p + geom_segment(aes(x=x,xend=xend,y=y,yend=yend), data=xf, size=1))
 } else {
   show(p)
 }
  ## pheatmap(ld,
  ##          cluster_rows=FALSE, cluster_cols=FALSE)
  invisible(ld)
}


bind.2.old <- function(x,y) {
  n <- nrow(x)
  m <- nrow(y)
  ni <- rep(1:n,each=m)
  mi <- rep(1:m,times=n)
  cBind(x[ni,,drop=FALSE], y[mi,,drop=FALSE])
}
l.rbind <- function(L, d1,d2) {
  do.call("rbind", lapply(L, function(x) x[[ d1 ]][[ d2 ]]))
}
l.c <- function(L,d1,d2) {
  do.call("c", lapply(L, function(x) x[[ d1 ]][[ d2 ]]))
}

myr2 <- function(X) {
  r2 <- ld(X,
           depth=ncol(X)-1,
           symmetric=TRUE,
           stat="R.squared")
  if(any(is.na(r2))) {
    r2.na <- as(is.na(r2),"matrix")
    use <- rowSums(r2.na)>0
    ## work around for r2=NA bug.  
    r2.cor <- as(cor(as(X[,use,drop=FALSE],"numeric"), use="pairwise.complete.obs")^2,"Matrix")
    r2[ which(r2.na) ] <- r2.cor[ which(r2.na[use,use]) ]
  }
  diag(r2) <- 1
  return(r2)
}

##' Derive tag SNPs for a SnpMatrix object using heirarchical clustering
##'
##' Uses complete linkage and the \code{\link{hclust}} function to define clusters, then cuts the tree at 1-tag.threshold
##' @title tag
##' @param snps colnames of the SnpMatrix object to be used
##' @param tag.threshold threshold to cut tree, default=0.99
##' @param samples optional, subset of samples to use
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @author Chris Wallace
tag <- function(X,tag.threshold=0.99, snps=NULL, samples=NULL) {
  if(!is.null(snps) || !is.null(samples))
    X <- X[samples,snps]
  r2 <- myr2(X)
   D <- as.dist(1-r2)
   hc <- hclust(D, method="complete")
   clusters <- cutree(hc, h=1-tag.threshold)
   snps.use <- names(clusters)[!duplicated(clusters)]
   r2.use <- r2[snps.use, colnames(X), drop=FALSE]
   tags <- rownames(r2.use)[apply(r2.use,2,which.max)]
   names(tags) <- colnames(r2.use)
   return(tags)
}

group.tags <- function(tags, keep) {
  groups <- tags[ names(tags) %in% tags.99 ]
  groups <- split(names(groups), groups)
}
##' Given a list of BMA results from glib, combine into a single object
##'
##' With the multicore package, it can be useful to farm out BMA jobs via mclapply
##' @title Collate a list of glib output
##' @param x a list of output from \code{\link{glib}}
##' @return an object similar to that returned by a single call of \code{glib}.  object$prior and object$posterior are _NOT_ correct
##' @author Chris Wallace
collate.bma <- function(results) {
  bma.mod <- results[[1]]
  bma.mod$inputs$models <- l.rbind(results, "inputs","models")
  bma.mod$models <- do.call("rbind", lapply(results, "[[", "models"))
  bma.mod$bf[[ "twologB10" ]] <- l.rbind(results, "bf","twologB10")
  bf.model <- exp(bma.mod$bf$twologB10/2)
  bma.mod$bf$postprob <-  t( t(bf.model) / colSums(bf.model))
  for(n in c("deviance","df","chi2","npar"))
    bma.mod$bf[[ n ]] <- l.c(results, "bf",n)
  for(n in c("mean","sd"))
    bma.mod$posterior.bymodel[[ n ]] <- l.c(results, "posterior.bymodel",n)
  for(n in c("coef","se"))
    bma.mod$glim.est[[ n ]] <- l.c(results, "glim.est",n)
  return(bma.mod)
}

##' Internal function: restrict to complete.cases of X and Y
##'
##' .. content for \details{} ..
##' @title make.data (internal function)
##' @param X matrix of explanatory variables
##' @param Y response
##' @param family character string, "binomial" or "gaussian", passed to glib
##' @return a list with entries X and Y
##' @author Chris Wallace
make.data <- function(X,Y,tags,family="binomial") {
  tags.names <- unique(tags)
  X <- X[,tags.names]
  if(is(X,"SnpMatrix") || is(X,"XSnpMatrix"))
    X <- as(X,"numeric")
  keep <- complete.cases(X)
  cat("Keeping",sum(keep),"of",nrow(X),"samples and",ncol(X),"SNPs\n")
  return(new("snpBMAdata",
             X=X[keep,],
             Y=Y[keep],
             tags=tags,
             family=family))
}



my.glib <- function(data, models) {
  if(data@family=="binomial") {
    error <- "binomial"
    link <- "logit"
  } else {
    error <- "gaussian"
    link <- "identity"
  }

  ## if(is.null(getOption("mc.cores"))) {
    cat("Evaluating",nrow(models),"models\n")
    bf <-   glib(data@X, data@Y, error=error, link=link, models=models, glimest=FALSE, post.bymodel=FALSE)$bf$twologB10
  ## } else {
  ##   cat("Evaluating",nrow(models),"models","over",getOption("mc.cores",1),"cores.\n")
  ##   njobs <- min(max(options()$mc.cores, 1), nrow(models)) # don't need more than 1 core per model
  ##   index <- rep(1:njobs, length=nrow(models))
  ##   inner.function <- function(j) {
  ##     glib(data@X, data@Y, error=error, link=link, models=models[j,,drop=FALSE], glimest=FALSE, post.bymodel=FALSE)$bf$twologB10
  ##   }
  ##   index.split <- split(1:nrow(models), index)
  ##   inner.results <- mclapply(index.split, inner.function)
  ##   ## collate results
  ##   bf <- do.call("rbind",inner.results)
  ##   models <- models[ unlist(index.split) , ]
  ## }  
  return(list(models=models,bf=bf))
}

bma.nsnps <- function(data, nsnps=1, groups=list(), models.drop=NULL) {
  if(nsnps<1 || nsnps>ncol(data@X))
    stop("nsnps must lie between 1 and ncol(data@X) inclusive")

  snps <- snps.use <- colnames(data@X)
  snps.exclude <- character(0)
  if(!is.null(models.drop)) {
    if(is.list(models.drop))
      models.drop <- stackModels(models.drop)
    snps.exclude <- snps(models.drop)
    snps.use <- setdiff(snps.use, snps.exclude)
    groups <- lapply(groups, setdiff, y=snps.exclude)
  }
  models <- make.models(snps.use, nsnps, groups=groups)
  nmodels <- nrow(models)
  if(!is.null(models.drop))
    models <- mdrop(models, models(models.drop))
  if(length(snps.exclude))
    models <- cbind2(models,
                    Matrix(0,nmodels,length(snps.exclude),
                           dimnames=list(NULL,snps.exclude),sparse=TRUE))[,snps]

  bma.run(data, models, nsnps, groups)

  ## x <- my.glib(data, models)
  ## return(new("snpBMA",
  ##            nsnps=nsnps,
  ##            nmodels=nmodels,
  ##            snps=data@tags,
  ##            groups=groups,
  ##            bf=x$bf,
  ##            models=x$models))             
}

bma.grow <- function(data, bma) {

  nsnps <- bma@nsnps + 1
  if(nsnps<1 || nsnps>ncol(data@X))
    stop("bma already contains maximum SNPs")
  models <- mgrow(bma)

  bma.run(data, models, nsnps, bma@groups)

}

bma.expand <- function(data, bma, groups) {

  nsnps <- bma@nsnps
  models <- mexpand(bma, groups)

  bma.run(data[,colnames(models)], models, nsnps, bma@groups)

}

## FIX!!
## > expand.snps<-c("rs2379078","rs4880781")
## > bma.e1 <- bma.expand(data.99, bma.1, groups=groups[expand.snps])
## groups not needed, creating a model matrix of 2 x 2 .
## Error in bind.2(models[i, -j], make.models.single(groups[[index.snp]],  : 
##   (list) object cannot be coerced to type 'double'

bma.run <- function(data, models, nsnps, groups) {

  x <- my.glib(data, models)
  
  return(new("snpBMA",
             nsnps=nsnps,
             nmodels=max.models(unique(data@tags), nsnps, groups),
             snps=data@tags,
             groups=groups,
             bf=x$bf,
             models=as(x$models,"dgCMatrix")))
}


## extend.matrix <- function(m, snps, by="row") {
##   dim.index <- if(by=="row") { 1 } else { 2 }
##   dimnames(m)[[ dim.index ]] <- unique(snps)
##   if(by=="row") {
##     m <- m[ snps , ]
##   } else {
##     m <- m[ , snps ]
##   }
##   dimnames(m)[[ dim.index ]] <- names(snps)
##   return(m)
## }

## extend.bytags <- function(object) {
##   object@models <- extend.matrix(object@models, object@snps, by="col")
##   return(object)
## }

drop.snps <- function(data, bma.result, bf.thr) {
  snps.drop <- unique(bma.result@snps)[bma.result@bf$twologB10[,2] < bf.thr]
  make.data(data@X[,setdiff(colnames(data@X),snps.drop)], data@Y,
            tags=setdiff(tags,snps.drop),family=data@family)
}
##' Binomial prior for number of SNPs in a model
##'
##' @title DELETE ME!
##' @param x number of SNPs in a model (defaults to 1:length(groups), ie returns a vector)
##' @param groups groups of SNPs, from which at most one SNP should be selected
##' @param expected expected number of SNPs in a model
##' @return prior probability/ies as a numeric
##' @export
##' @author Chris Wallace
prior.nsnps <- function(x=1:length(groups), groups, expected) {
  n <- length(groups)
  p <- expected/n
  dbinom(x, size=n, prob=p)  
}

bma.combine <- function(L) {
  ## L <- lapply(L, extend.bytags)
  models <- lapply(L, function(l) l@models)
  bf <- c(list(1),lapply(L, function(l) l@bf$twologB10))
  names(bf) <- 0:length(L)
  bf.sum <- sapply(bf,sum)
  print(posterior.nsnps <- bf.sum/sum(bf.sum))
  coef <- lapply(L, function(l) do.call("rbind", l@glim.est$coef))
  se <- lapply(L, function(l) do.call("rbind", l@glim.est$se)) 

  ## df1 <- data.frame(snp=unique(L[[1]]@snps),
  ##                   bf.1=bf[[1]][,2])
  ## index <- which(models[[2]]==1, arr.ind=TRUE)
  ## index <- t(apply(models[[2]]==1, 1, which))
  ## df2 <- as.data.frame(index)
  ## df2 <- cbind(df2, lapply(df2, function(i) unique(L[[2]]@snps)[i]))
  ## colnames(df2) <- c("i1","i2","snp1","snp2")
  ## df2$bf1 <- df2$bf2 <- bf[[2]][,2]
  ## df <- merge(df1,df2[,c("snp1","bf1")], by.x="snp", by.y="snp1")
  ## df <- merge(df,df2[,c("snp2","bf2")], by.x="snp", by.y="snp2")
  ## ggplot(df, aes(x=bf1,y=bf2,col=bf.1<0)) + geom_point()
  
  ## library(ggplot2)
  ## ggplot(df,aes(x=bf.1,y=bf.2)) + geom_point()

}

################################################################################

## old stuff from groups


##' create groups of SNPs with LD > r2.threshold with index.snps
##'
##' .. content for \details{} ..
##' @title index.groups
##' @param X SnpMatrix object
##' @param index.snps snps to index groups
##' @param r2.threshold threshold above which snps should be grouped with an index snp
##' @param snps SNP subset to be used (columns of X)
##' @param samples sample subset to be used (rows of X)
##' @return list of character vectors, each vector being a single group and the names of the list being the index snps
##' @author Chris Wallace
index.groups <- function(X, index.snps, r2.threshold=0.9, snps=NULL, samples=NULL) {
  if(r2.threshold<0 || r2.threshold>1)
    stop("r2.threshold must lie within (0,1)")
  index.snps <- unique(index.snps)
  if(!is.null(snps))
    snps <- unique(snps)
  if(!is.null(snps) & !is.null(samples)) {
    X <- X[ samples, snps ]
  } else {
    if(!is.null(snps))
      X <- X[,snps]
    if(!is.null(samples))
      X <- X[samples,]
  }
  ld.groups <- as.data.frame(t(ld(X[,index.snps],
                                  X,
                                  stat="R.squared")))
  group.snps <- lapply(ld.groups, function(x)
                       rownames(ld.groups)[which(x>=r2.threshold)])
  ## ## merge groups in LD
  ## if("imm_10_6137289" %in% names(group.snps))
  ##   group.snps <- do.merge(group.snps, c("imm_10_6137289","imm_10_6134703"))

  ## if("imm_10_6142018" %in% names(group.snps))
  ##   group.snps <- do.merge(group.snps, c("imm_10_6142018","ccc.10.6103325.C.T"))

  ## if("imm_10_6116082" %in% names(group.snps))
  ##   group.snps <- do.merge(group.snps, c("imm_10_6116082","imm_10_6081860"))

  return(group.snps)
}

##' create groups of SNPs with LD > r2.threshold with index.snps
r2.groups <- function(X,r2.threshold=0.9, snps=NULL, samples=NULL, do.plot=FALSE) {
  if(!is.null(snps) || !is.null(samples))
    X <- X[samples,snps]
  r2 <- myr2(X)
  D <- as.dist(1-r2)
  hc <- hclust(D, method="complete")
  clusters <- cutree(hc, h=1-r2.threshold)
  if(do.plot) {
    plot(hc)
    rect.hclust(hc, h=1-r2.threshold)
  }
  groups <- split(names(clusters),clusters)
  names(groups) <- sapply(groups,"[[",1)
  return(groups)
}



## do.merge <- function(snps, merge.pair) {
##   if(!all(merge.pair %in% names(snps)))
##     stop("merge.pair must index snps\n")
##   merged.group <- unique(unlist(snps[merge.pair]))
##   snps[[ merge.pair[1] ]] <- merged.group
##   wh.drop <- which(names(snps)==merge.pair[2])
##   return(snps[-wh.drop])
## }

## ##' Internal function: restrict to complete.cases of X and Y
## ##'
## ##' .. content for \details{} ..
## ##' @title make.data (internal function)
## ##' @param X matrix of explanatory variables
## ##' @param Y response
## ##' @return a list with entries X and Y
## ##' @author Chris Wallace
## make.data <- function(X,Y) {
##   keep <- complete.cases(X)
##   cat("Keeping",sum(keep),"of",nrow(X),"samples and",ncol(X),"SNPs\n")
##   return(list(X=X[keep,], Y=Y[keep]))
## }

## make.models <- function(n.group) {
##  models.list <- lapply(n.group, function(n) {
##     mat <- matrix(0,nrow=n+1,ncol=n)
##     diag(mat) <- 1
##     mat
##   })
##   models <- models.list[[1]]
##   if(length(models.list)>1) {
##     for(i in 2:length(models.list))
##       models <- bind.2(models,models.list[[i]])
##   }
##  return(models)
## }

## drop.groups <- function(group.result, to.drop=group.result$snps) {

##   stop("DOESN'T WORK YET!!!\n")
  
##   if(length(to.drop)==0)
##     stop("This function assesses evidence for dropping (groups of) of SNPs.  to.drop must have length>=1.\n")
##   if(!all(unlist(to.drop) %in% unlist(group.result$snps)))
##     stop("All to.drop must be listed in snps.\n")
  
##   ## bayes factors for merging or not
##   ## merge if BF(both) << BF(either)
##   models <- group.result$bma.mod$models
##   colnames(models) <- unlist(group.result$snps)
##   bf.in.out <- matrix(NA,length(to.drop),3)
##   rownames(bf.in.out) <- sapply(to.drop,paste,collapse="-")
##   tbf.mod <- t(exp(group.result$bma.mod$bf$twologB10[,2]/2))
##   for(i in 1:length(to.drop)) {
##     idrop <- to.drop[[i]]
##     mod.in <- which(apply(models[, idrop, drop=FALSE ]==1, 1, any))
##     mod.out <- which(apply(models[, idrop, drop=FALSE ]==0, 1, all))
##     bf.in <- 2 * log( tbf.mod[,mod.in] %*% models[mod.in, ])
##     bf.out <- 2 * log( tbf.mod[,mod.out] %*% models[mod.out, ])
##     pp.in <- exp(bf.in/2) / sum(exp(bf.in/2))
##     pp.out <- exp(bf.out/2) / sum(exp(bf.out/2))
##     pp.out <- 2 * log( tbf.mod[,mod.out] %*% models[mod.out, ])
##     bf.out <- group.result$bma.mod$bf$twologB10[ mod.out, ]
##     bf.in.out[i,] <- 2*(log(colSums( exp(bf.in/2) )) - log(colSums( exp(bf.out/2) )))
##   }
  
##   cat("twologBF for including SNPs vs not. Values > 0 favour including both SNPs in the model.\n")
##   print(bf.in.out)
  
## }

## check.merge <- function(group.result, merge.candidates) {

##   if(length(merge.candidates)!=2)
##     stop("This function assesses evidence for merging _pairs_ of SNPs.  merge.candidates must have length==2.\n")
##   if(!all(merge.candidates %in% unlist(group.result$snps)))
##     stop("All merge.candidates must be listed in snps.\n")
  
##   ## bayes factors for merging or not
##   ## merge if BF(both) << BF(either)
##   models <- group.result$bma.mod$models
##   colnames(models) <- unlist(group.result$snps)
##   mod.both <- which(models[, merge.candidates[1] ]==1 &
##                     models[, merge.candidates[2] ]==1)
##   mod.either <- setdiff(which(models[, merge.candidates[1] ]==1 |
##                               models[, merge.candidates[2] ]==1),
##                         mod.both)
##   bf.both <- group.result$bma.mod$bf$twologB10[ mod.both, ]
##   bf.either <- group.result$bma.mod$bf$twologB10[ mod.either, ]
##   bf.both.either = 2*(log(colSums( exp(bf.both/2) )) - log(colSums( exp(bf.either/2) )))
##   cat("twologBF for including both SNPs vs either. Values > 0 favour including both SNPs in the model.\n")
##   print(bf.both.either)
##   invisible(bf.both.either)
 
## }

## ##' BMA for SNPs in a single group with bootstrapping
## ##'
## ##' Allow exactly one SNP from the group to be included
## ##' @title boot.bma
## ##' @param snps character vector of SNPs describing the model space to be explored
## ##' @param tag.threshold SNPs with r2>tag.threshold will be grouped, and just one used to tag the others
## ##' @param B number of bootstraps.  Set B=0 for no bootstraps.
## ##' @param additional.snps a character vector of additional SNPs which will be included in all models, but see also \code{group.results}
## ##' @return a list with three elements: BMA results, snps used and tags
## ##' @author Chris Wallace
## boot.bma <- function(snps, tag.threshold=0.999, B=100, additional.snps=NULL) {
##   ## consider all pairs of (linked) snps in a group.  What is the BF between every pair?  What is its associated uncertainty?
##   cat(snps,"\n")
##   if(length(snps)==1)
##     return(NULL)
##   all.snps <- c(snps,additional.snps)
##   X.sub <- X[,all.snps]
##   keep <- complete.cases(X.sub)
##   cat("Keeping",sum(keep),"of",nrow(X.sub),"samples and",length(snps),"SNPs\n")
##   X.sub <- X.sub[keep,]
##   Y.sub <- Y[keep]
##   ## ld thin by minimum possible
##   tags <- tag(X=snp.data, snps=snps, tag.threshold=tag.threshold)
##   snps.use <- unique(tags)
##   n <- length(snps.use)
##   models <- matrix(0,n,n)
##   models[ cbind(1:n, 1:n) ] <- 1
##   if(nrow(wh)) {
##     for(i in 1:nrow(wh)) {
##       models <- models[ !(models[,1]==wh[i,1] & models[,2]==wh[i,2]), ]
##     }
##   }

##   models <- cbind(models, matrix(1,nrow(models), length(additional.snps)))      
##   data <- cbind(Y.sub,X.sub[,c(snps.use, additional.snps)])

##   boot.fn <- function(data,i=1:nrow(data)) {
##     ystar <- data[i,1]
##     xstar <- data[i,-1]
##     bma.star <- try(glib(xstar, ystar, error="binomial", link = "logit",
##                      models=models), silent=TRUE)
##     if(inherits(bma.star, "try-error"))
##       return(rep(NA,n))
##     return(bma.star$bf$twologB10[1:n,2])
##   }  
  
##   bma.boot <- boot(data, boot.fn, B)
##   return(list(bma.boot=bma.boot, tags=tags, snps=snps))

## }
## ##' BMA for SNPs in a single group, priors determined from a previous BMA
## ##'
## ##' Allow exactly one SNP from the group to be included
## ##' @title within.bma
## ##' @param snps character vector of SNPs describing the model space to be explored
## ##' @param tag.threshold SNPs with r2>tag.threshold will be grouped, and just one used to tag the others
## ##' @param additional.snps a character vector of additional SNPs which will be included in all models, but see also \code{group.results}
## ##' @param group.results if group.results is given (the output of \code{group.bma()}), different combinations of additional.snps will be considered and the prior of each group determined from group.bma
## ##' @return a list with 4 elements: bma.results, priors, tags, snps
## ##' @author Chris Wallace
## within.bma <- function(snps, tag.threshold=0.99, additional.snps=NULL,
##                        group.results=NULL, family="binomial",
##                        X, Y) {
##   ## consider all pairs of (linked) snps in a group.  What is the BF between every pair?  What is its associated uncertainty?
##   cat(snps,"\n")
##   if(length(snps)==1)
##     return(NULL)
##   if(!is.null(additional.snps))
##     additional.snps <- setdiff(additional.snps, snps)
##   all.snps <- c(snps,additional.snps)
##   data <- make.data(X[,all.snps], Y)
  
##   ## ld thin by minimum possible
##   tags <- tag(snps=snps, tag.threshold=tag.threshold, samples=rownames(data$X))
##   snps.use <- unique(tags)
##   n <- length(snps.use)
##   cat("Using",n,"SNPs to tag",length(snps),"at r2 >=",tag.threshold,"\n")
##   models <- matrix(0,n,n)
##   diag(models) <- 1
  
##   priors <- 1
##   models2 <- NULL
  
##   if(!is.null(additional.snps) & is.null(group.results)) {
##     models2 <- matrix(1,1, length(additional.snps))
##   }
##   if(!is.null(additional.snps) & !is.null(group.results)) {
##     index.snps <- names(group.results$snps)
##     test.index <- which(index.snps %in% snps)
##     if(length(test.index)!=1)
##       stop(paste("Exactly one test snp must index a group in group.results.  Found",
##                  length(test.index),"\n"))
##     group.models <- do.call("rbind",
##                             lapply(strsplit(rownames(group.results$gmod.result),""),
##                                    as.integer))    
##     use <- which(group.models[,test.index]==1)
##     priors <- group.results$mod.result[use,"postprob"]
##     priors <- priors/sum(priors)
##     models2 <- group.models[use,-test.index]
##   }

##   data$X <- data$X[,c(snps.use, additional.snps)]

##   if(is.null(group.results) || is.null(getOption("mc.cores"))) {
##     cat("Evaluating",nrow(models),"models\n")
##     bma.results <- my.glib(data,  family=family,
##                      models=models)
##   } else {
##     njobs <- min(max(options()$mc.cores, 1), nrow(models2)) # don't need more than 1 core per model
##     index <- rep(1:njobs, length=nrow(models2))
##     inner.function <- function(j) {
##       wh <- which(index==j)
##       inner.results <- lapply(as.list(wh), function(k) {
##         tmp <- my.glib(data,  family=family,
##                        models=cbind(models,
##                          matrix(models2[k,],nrow=nrow(models),ncol=ncol(models2),byrow=TRUE)))
##         tmp$bf$twologB10 <- tmp$bf$twologB10 + 2 * log( priors[k] )
##         return(tmp)
##       })
##       ## collate results
##       collate.bma(inner.results)
##     }
##     cat("Evaluating",nrow(models),"models x",length(priors),"priors ( total =",nrow(models)*length(priors),") over",getOption("mc.cores"),"cores.\n")
##     bma.results <- collate.bma(mclapply(1:njobs, inner.function))
##     priors <- unlist(lapply(1:njobs, function(j) priors[ which(index==j) ]))
##   }
  
##   return(list(bma.results=bma.results, priors=priors, tags=tags, snps=snps))

## }


## ##' Summarize a boot.bma result
## ##'
## ##' .. content for \details{} ..
## ##' @title summ.boot
## ##' @param x object returned by \code{\link{boot.bma}}
## ##' @return data.frame
## ##' @author Chris Wallace
## summ.boot <- function(x) {
##   if(is.null(x))
##     return(NULL)
##   b <- x$bma.boot
##   snp.best <- which.max(b$t0)
##   this.r2 <- r2[x$snps[snp.best],x$snps]
##   m <- match(x$tags,x$snps)
##   bf.diff <- (b$t0 - b$t0[snp.best])[m]
##   bf.diff.star <- (b$t - b$t[,snp.best])[,m]
##   n <- length(bf.diff)
##   summ <- data.frame(snp=x$snps,
##                      r2=this.r2,
##                      bf.diff,
##                      lq=apply(bf.diff.star,2,quantile,probs=0.25, na.rm=TRUE),
##                      uq=apply(bf.diff.star,2,quantile,probs=0.75, na.rm=TRUE),
##                      lci=apply(bf.diff.star,2,quantile,probs=0.025, na.rm=TRUE),
##                      uci=apply(bf.diff.star,2,quantile,probs=0.975, na.rm=TRUE))
## }

## summ.within <- function(x) {
##   if(is.null(x))
##     return(NULL)
##   b <- x$bma.results 
##   n <- length(unique(x$tags))
##   pr <- x$priors
##   bf.full <- matrix(b$bf$twologB10[,2],ncol=n,byrow=TRUE)
##   colnames(bf.full) <- unique(x$tags)
##   bf.full <- bf.full[, x$tags]
##   colnames(bf.full) <- names(x$tags)
##   ## bf.full <- t(sapply(b, function(x) x$bf$twologB10[,2]))
##   bf.marg <- as.vector(t(pr) %*% exp(bf.full/2))
##   models <- b$models[,1:n]
##   pp <- bf.marg/sum(bf.marg)
##   bf.marg <- 2 * log(bf.marg)
##   snp.best <- which.max(bf.marg)
##   ctrl.r2 <- as.vector(ld(snp.data[sample.support$cc==0,x$snps[snp.best],drop=FALSE],
##                 snp.data[sample.support$cc==0,x$snps],
##                 stat=c("R.squared"),symmetric=TRUE))
##   case.r2 <- as.vector(ld(snp.data[sample.support$cc==1,x$snps[snp.best],drop=FALSE],
##                 snp.data[sample.support$cc==1,x$snps],
##                 stat=c("R.squared"),symmetric=TRUE))
  
##   bf.diff <- bf.marg - bf.marg[snp.best]
##   summ <- data.frame(snp=x$snps,
##                      r2=ctrl.r2,
##                      r2.cse=case.r2,
##                      bf.marg,
##                      pp,
##                      bf.diff)
##   return(summ)
## }

## best.within <- function(x,bf.threshold=-6) {
##   summ <- lapply(x, summ.within)
##   use <- !sapply(summ,is.null)
##   for(i in 1:length(summ))
##     summ[[i]]$group <- i
##   summ <- do.call("rbind",summ[use])   
##   subset(summ, bf.diff > -6)
## }

## plot.within <- function(x) {
##    summ <- lapply(x, summ.within)
##    use <- !sapply(summ,is.null)
##    for(i in 1:length(summ))
##      summ[[i]]$group <- i
##    summ <- do.call("rbind",summ[use])   
##    ggplot(summ, aes(x=r2,y=bf.diff,col=as.factor(group))) +
##     geom_hline(yintercept=c(0,-2.2,-6,-10),col="grey") +
##       geom_point() +
##         geom_line() +
## #    geom_text(x=1,y=0,label=x$snps[snp.best],angle=315,hjust=1,vjust=1) +
##         geom_text(aes(label=snp),angle=315,hjust=1,vjust=1) +
##         xlab("R.squared") + ylab("2logBF")
## }
 

## plot.boot <- function(x) {
##   summ <- summ.boot(x)
##    ggplot(summ, aes(x=r2,y=bf.diff,ymin=lq,ymax=uq)) +
##     geom_hline(yintercept=c(0,-2.2,-6,-10),col="grey") +
##       geom_pointrange() +
## #    geom_text(x=1,y=0,label=x$snps[snp.best],angle=315,hjust=1,vjust=1) +
##         geom_text(aes(label=snp),angle=315,hjust=1,vjust=1) +
##         xlab("R.squared") + ylab("2logBF")
## }

## ##' BMA for SNPs in groups
## ##'
## ##' Allow at most one SNP from each group in the model
## ##' @title group.bma
## ##' @param group.snps list of character vectors, each vector describes the SNPs in one group
## ##' @return list of results=bma.results,
## ##' \itemize{
## ##'  \item{results}{snp and group twologB10 and posterior probs}
## ##'  \item{gmod.result}{group level model results}
## ##'  \item{mod.result}{full model results}
## ##' \item{effects}{mean, sd of effects in each model. models ordered as in mod.result}
## ##' }
## ##' @author Chris Wallace
## group.bma <- function(group.snps, tag.threshold=0.99, family="binomial", X=X, Y=Y) {
##   ## average over all snps in a group - what is rel PIP of a group?
##   if(!is.list(group.snps)) ## assume want groups with a single snp in each
##     group.snps <- structure(as.list(group.snps), names=group.snps)
##   if(!all(unlist(group.snps) %in% colnames(X)))
##     stop("Some group.snps not found in X.\n")
##   if(is.null(names(group.snps)) & all(sapply(group.snps,length)==1))
##     names(group.snps) <- unlist(group.snps)
##   data <- make.data(X[,unlist(group.snps)], Y)

##   ## ld thin by minimum possible
##   tags <- tag(X=snp.data, snps=unlist(group.snps), tag.threshold=tag.threshold, samples=rownames(data$X))
##   group.use <- lapply(group.snps, intersect, tags)
##   n.group <- sapply(group.use,length)
##   n.total <- sum(n.group)
##   group.num <- split(1:n.total, rep(1:length(group.use),times=n.group))
##   cat("Keeping",n.total,"SNPs of",sum(sapply(group.snps,length)),"after thinning at r2=",tag.threshold,"\n")
##   data$X <- data$X[, unlist(group.use)]
  
##   models <- make.models(n.group)
    
##   if(!is.null(getOption("mc.cores"))) {
    
##     njobs <- min(max(options()$mc.cores, 1), nrow(models)) # don't need more than 1 core per model
##     cat("Evaluating",nrow(models),"models across",njobs,"cores.\n")
##     index <- rep(1:njobs, length=nrow(models))
##     inner.function <- function(j) {
##       my.glib(data=data, family=family,
##            models=models[index==j,,drop=FALSE])
##     }
##     results <- mclapply(1:njobs, inner.function)
    
##     ## collate results
##     bma.mod <- collate.bma(results)

##   } else {  

##     cat("Evaluating",nrow(models),"models on a single core.\n")
##     bma.mod <- my.glib(data=data, family=family,
##                     models=models)

##   }
  
##   ## collapse models by group
##   gmodels <- sapply(group.num, function(g) {
##     rowSums(bma.mod$inputs$models[, g, drop=FALSE])
##   })
##   bf.model <- bma.mod$bf$twologB10[,2]
##   bf.group <- 2*log(t(exp(bf.model/2)) %*% gmodels)
##   pp.model <- bma.mod$bf$postprob[,2]
##   pp.group <- t(pp.model) %*% gmodels
##   pp.model <- tapply(pp.model, apply(gmodels,1,paste,collapse=""), sum)
##   o<-order(pp.model, decreasing=TRUE)
##   group.results <- data.frame(group=rep(1:length(n.group)),
##                               pp.group=t(pp.group),
##                               bf.group=t(bf.group))
##   gmod.results <- data.frame(postprob=pp.model[o], twologB10=bf.model[o])

##   ## summarize by snp
##   pp.snp <- t(bma.mod$bf$postprob[,2]) %*% bma.mod$inputs$models
##   bf.snp <- 2 * log(t(exp(bma.mod$bf$twologB10[,2]/2)) %*% bma.mod$inputs$models)
##   snp.results <- data.frame(snp=unlist(group.use),
##                             group=rep(1:length(n.group), times=n.group),
##                             pp.snp=t(pp.snp),
##                             bf.snp=t(bf.snp))
##   snp.results$tag <- tags[ unlist(group.use) ]
  
##   mod.results <- data.frame(postprob=bma.mod$bf$postprob[,2],
##                             twologB10=bma.mod$bf$twologB10[,2])
##   rownames(mod.results) <- apply(bma.mod$inputs$models,1,paste,collapse="")

##   bma.results <- merge(snp.results,group.results)
 
##   list(results=bma.results,
##        gmod.results=gmod.results,
##        mod.result=mod.results,
##        effects=bma.mod$posterior.bymodel,
##        snps=group.snps,
##        tags=tags,
##        bma.mod=bma.mod)       
## }
