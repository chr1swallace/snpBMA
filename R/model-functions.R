make.models.slow <- function(snps,n.use,groups=NULL) {
  if(is.null(groups) || n.use==1) { # simple
    combs <- combn(length(snps),n.use)
    models <- matrix(0,ncol(combs),length(snps))
    colnames(models) <- snps
    ind <- cbind(as.vector(col(combs)),as.vector(combs))
    models[ind] <- 1
  } else {
    n.groups <- sapply(groups,length)
    if(length(snps) != sum(n.groups) || !all(snps %in% unlist(groups)))
      stop("all snps should be in a single group")
    models.list <- lapply(n.groups, function(n) {
      diag(nrow=n,ncol=n)
    })
    null.list <- lapply(n.groups, function(n) {
      matrix(0,nrow=1,ncol=n)
    })
    ## pick n.use from models.list, and remainder from null.list, combine
    combs <- combn(length(n.groups),n.use)
    models <- vector("list",ncol(combs))
    for(j in 1:ncol(combs)) {
      models[[j]] <- outer.models(models.list[[combs[1,j]]], models.list[[combs[2,j]]])
      if(nrow(combs)>2) {
        for(i in 3:nrow(combs))
          models[[j]] <- outer.models(models[[j]], models.list[[ combs[i,j] ]])
      }
      models[[j]] <- cbind(models[[j]], matrix(0,nrow(models[[j]]), sum(n.groups[ -combs[,j] ])))
      colnames(models[[j]]) <- c(unlist(groups[ combs[,j] ]), unlist(groups[ -combs[,j] ]))
      models[[j]] <- models[[j]][,snps]
    }
    models <- do.call("rbind",models)
  }
  return(models)
}

outer.models <- function(x,y=NULL) {
  if(!is.null(y))
    x <- list(x,y)

  snps <- unlist(lapply(x, colnames))
  nr <- lapply(x,nrow)
  if(any(nr==0)) # deal with null models
    return(null.model(snps=snps))

  ind <- expand.grid(lapply(nr, function(n) if(n>0) { 1:n } else { 0 }))
  
  nc <- sapply(x,ncol)
  beg <- c(1,cumsum(nc)[-length(nc)]+1)
  end <- cumsum(nc)
  if(all(sapply(x,"class")=="matrix")) {
    mat <- matrix(0,nrow=nrow(ind),ncol=sum(nc))
  } else {
    mat <- Matrix(0,nrow=nrow(ind),ncol=sum(nc),sparse=TRUE)
  }
  for(i in 1:length(ind))
    if(nrow(x[[i]])) # deal with null models
      mat[, beg[[i]]:end[[i]] ] <- x[[i]][ ind[[i]], ]
##  return(Matrix(mat,sparse=TRUE))

  colnames(mat) <- snps
  return(mat)
##  system.time( do.call("cBind",mapply(function(mat,i) { mat[i,,drop=FALSE] }, x, as.list(ind), SIMPLIFY=FALSE)))
}
##' Internal function, null model matrix
##'
##' @title null.model
##' @param snps character vector of snp names
##' @return Matrix, nrow=0, ncol=length(snps), colnames=snps
##' @author Chris Wallace
null.model <- function(snps=NULL,nsnps=length(snps)) {
  M <- Matrix(0,0,nsnps)
  if(!is.null(snps))
    colnames(M) <- snps
  return(M)
}
##' Internal function, return a model matrix for choosing n.use snps out of snps
##'
##' @title make.models.single
##' @param snps character vector of snp names
##' @param n.use integer: number of snps to use in each model
##' @param quiet if TRUE, suppress progress messages
##' @return Matrix, with nrow=number of models, ncol=number of SNPs,
##' entries of 1/0 indicate whether SNP is in/excluded from model
##' @author Chris Wallace
make.models.single <- function(snps,n.use,quiet=FALSE) {
  if(n.use > length(snps))
    return(null.model(snps))

  combs <- combn(length(snps),n.use)
  if(!quiet)
    cat("groups not needed, creating a model matrix of",ncol(combs),"x",length(snps),".\n")
  models <- Matrix(0,ncol(combs),length(snps),sparse=TRUE)
  colnames(models) <- snps
  ind <- cbind(as.vector(col(combs)),as.vector(combs))
  models[ind] <- 1
  return(models)
}
##' Internal function, return a model matrix for choosing n.use snps out of snps
##'
##' @title make.models.multi
##' @param snps character vector of snp names
##' @param n.use integer: number of snps to use in each model
##' @param groups list of character vectors, each vector indicating a
##' group of SNPs from which a maximum of one may be selected
##' @param quiet if TRUE, suppress progress messages
##' @return  Matrix, with nrow=number of models, ncol=number of SNPs,
##' entries of 1/0 indicate whether SNP is in/excluded from model
##' @author Chris Wallace
make.models.multi <- function(snps,n.use,groups,quiet=FALSE) {
  if(length(groups) < n.use)
    return(null.model(snps))
  
  n.groups <- sapply(groups,length)
  
  combs <- combn(length(n.groups),n.use)
  models.list <- lapply(n.groups, function(n) {
      diag(nrow=n,ncol=n)
  })
  models <- vector("list",ncol(combs))
  for(j in 1:ncol(combs)) {
    ##        cat(".")
    models[[j]] <- outer.models(models.list[ combs[,j] ])
    models[[j]] <- cbind(models[[j]],
                         matrix(0,nrow(models[[j]]), sum(n.groups[ -combs[,j] ])))
    colnames(models[[j]]) <- c(unlist(groups[ combs[,j] ]), unlist(groups[ -combs[,j] ]))
    models[[j]] <- models[[j]][,snps]
  }
  nr <- sapply(models,nrow)
  models.final <- Matrix(0,sum(nr),ncol(models[[1]]))
  colnames(models.final) <- snps
  end <- cumsum(nr)
  beg <- c(1,end[-length(end)]+1)
  for(j in 1:length(models)) {
    models.final[beg[[j]]:end[[j]],] <- models[[j]]
  }
  ##      models <- do.call("rBind",models)
  return(models.final)
}

model.empty <- function(snps) {
  Matrix(0,nrow=1,ncol=length(snps),sparse=TRUE,
         dimnames=list(NULL,snps))
}

make.models <- function(snps,n.use,groups=list(),quiet=FALSE) {
  if(!length(groups) || n.use==1) # simple
    return(make.models.single(snps,n.use,quiet))
  
  n.groups <- sapply(groups,length)
  if(length(snps) != sum(n.groups) || !all(snps %in% unlist(groups)))
    stop("all snps should be in a single group")
  
  if(all(n.groups==1)) ## all singletons, ignore grouping
    return(make.models.single(snps,n.use,quiet=TRUE))
  
  if(all(n.groups>1)) ## all multi, generate models using at most one SNP from each group
    return(make.models.multi(snps,n.use,groups,quiet=TRUE))
  
  ## mixed single and multi
  ## deal with single snp groups first - these are easy
  which.single <- which(n.groups==1)
  if(!quiet)
    cat("mixture of",length(groups)-length(which.single),"groups and",length(which.single),"singles.\n")
  if(!quiet)
    cat("singles first...\n")
  single.list <- c(list(model.empty(unlist(groups[which.single]))),
                   lapply(as.list(1:n.use), function(i) {
                     make.models(unlist(groups[which.single]), n.use=i,quiet=TRUE) }))
  ## now each of these needs to outer.models with all n:0 possibilities from remainder
  if(!quiet)
    cat("multis next...\n")
  multi.list <- c(lapply(as.list(n.use:1), function(i) {
    make.models(unlist(groups[-which.single]), n.use=i, groups=groups[-which.single], quiet=TRUE)
  }),
                  list(model.empty(unlist(groups[-which.single]))))
  if(!quiet)
    cat("outer.models\n")
  models.list <- mapply(outer.models, single.list, multi.list) ## mcmapply
  if(!quiet)
    cat("finally rBind\n")
  models <- do.call("rBind",models.list)
  colnames(models) <- c(unlist(groups[which.single]),unlist(groups[-which.single]))
  models <- models[,snps]
  return(models)
}

##' Compare potential parent and child models according to BF
##'
##' Returns models whose children SHOULD NOT be visited
##' @title mcomp
##' @param parents snpBMA object
##' @param children snpBMA object
##' @param ... arguments passed to mcomp.detail()
##' @return a Matrix containing all child models for which twologBF (parent/child) > 2*log(rel.prior)
##' @author Chris Wallace
models.diff <- function(parents, children, ...) {
  ## for each child model, identify its parents
  ## models to drop should be defined as the set with any 2*logbf(parent/child) > 2*log(rel.prior) + 2*lbf
  return(mcomp.detail(m.parent=parents@models, m.child=children@models,
                      bf.parent=parents@bf, bf.child=children@bf,
                      ntarget=parents@nsnps, what="drop", ...))
}

##' Compare potential parent and child models according to BF
##'
##' Returns models whose children SHOULD be visited
##' @title mcomp
##' @param bma snpBMA object
##' @param snps character vector of snps to drop
##' @param ... arguments passed to mcomp.detail()
##' @return a snpBMA object with a subset of models in the input
##' \code{children} object, all child models for which twologBF
##' (parent/child) > 2*log(pp.fold)
##' @author Chris Wallace
models.prune.snps <- function(bma, snps.to.drop, ...) {
  ## for each child model, identify its parents
  ## models to drop should be defined as the set with any 2*logbf(parent/child) > 2*log(rel.prior) + 2*lbf
  if(!all(snps.to.drop %in% snps(bma)))
    stop("SNPs to drop must be in BMA object\n")

  models <- bma@models
  models.drop <- which(rowSums(models[,snps.to.drop]) > 0)
  if(length(models.drop))
    return(bma[-models.drop, ])

  return(bma)
}
 
##' Compare potential parent and child models according to BF
##'
##' Returns models whose children SHOULD be visited
##' @title mcomp
##' @param bma snpBMA object
##' @param snps character vector of snps to drop
##' @param ... arguments passed to mcomp.detail()
##' @return a snpBMA object with a subset of models in the input
##' \code{children} object, all child models for which twologBF
##' (parent/child) > 2*log(pp.fold)
##' @author Chris Wallace
snps.prune.data <- function(data, snps.to.drop, ...) {
  ## for each child model, identify its parents
  ## models to drop should be defined as the set with any 2*logbf(parent/child) > 2*log(rel.prior) + 2*lbf
  tags <- data@tags
  tags <- tags[!(tags %in% snps.drop)]

  make.data(data@X, data@Y,tags=tags, family=data@family)
}

##' Compare potential parent and child models according to BF
##'
##' Returns models whose children SHOULD be visited
##' @title mcomp
##' @param parents snpBMA object
##' @param children snpBMA object
##' @param ... arguments passed to mcomp.detail()
##' @return a snpBMA object with a subset of models in the input
##' \code{children} object, all child models for which twologBF
##' (parent/child) > 2*log(pp.fold)
##' @author Chris Wallace
models.prune <- function(parents, children, ...) {
  ## for each child model, identify its parents
  ## models to drop should be defined as the set with any 2*logbf(parent/child) > 2*log(rel.prior) + 2*lbf
  newmodels <- mcomp.detail(m.parent=parents@models, m.child=children@models,
                            bf.parent=parents@bf, bf.child=children@bf,
                            ntarget=parents@nsnps, n.child=children@nsnps, what="keep", ...)
  return(bma[ newmodels, ])
}
##' mcomp.detail, internal function
##'
##' @title mcomp.detail
##' @param m.parent 
##' @param m.child 
##' @param bf.parent 
##' @param bf.child 
##' @param ntarget 
##' @param n.child 
##' @param prior.parents 
##' @param prior.children 
##' @param what 
##' @param pp.fold 
##' @param quiet 
##' @return object of class dropModels, or an index vector of which
##' rows of supplied m.child should be dropped
##' @author Chris Wallace
mcomp.detail <- function(m.parent, m.child, bf.parent, bf.child, ntarget, n.child, prior.parents, prior.children, what, pp.fold=10, quiet=FALSE) {
  ## for each child model, identify its parents
  ## models to drop should be defined as the set with any 2*logbf(parent/child) > 2*log(rel.prior) + 2*lbf

  cols.use <- intersect(colnames(m.parent),colnames(m.child))
  relate <- m.parent[,cols.use] %*% Matrix::t(m.child[,cols.use])
##  print(relate[1:10,1:10])
  index <- Matrix::which(relate==ntarget, arr.ind=TRUE)
  bf.parent <- bf.parent[ index[,1], 2 ]
  bf.child <- bf.child[ index[,2], 2 ]
  relate2 <- Matrix(0,nrow(relate),ncol(relate),sparse=TRUE)
  relate2[index] <- bf.parent + 2*log(prior.parents) - bf.child - 2*log(prior.children) - 2*log(pp.fold)
  ## rel.prior is prior P(M|child)/P(M|parent)
  cmp2 <- relate2 != 0
  cmp1 <- relate2 > 0
  ##child.all <- apply(cmp1 | !cmp2, 2, all)
  child.any <- apply(cmp1 & cmp2, 2, any)
  if(!quiet)
    cat("Identified",sum(child.any),"of",nrow(m.child),
        "models with pp(parent) >", pp.fold,"* pp(child)\n")

  if(what=="keep") { ## return an index of models to keep
    return(!child.any)
  }
  
  models.keep <- new("Models",models=m.child[!child.any,,drop=FALSE], nsnps=n.child)  
  if(what=="drop") { ## return a matrix of models to drop
    models.drop <- m.child[child.any,,drop=FALSE ]
    snps.drop <- apply(models.keep==0,2,all)
    return( new("dropModels", models=models.drop, snps=colnames(m.child)[snps.drop] ) )
  }
}

mdrop <- function(models, drop, quiet=FALSE) {
  if(!quiet) {
    cat("Dropping models:\n")
    show(drop)
  }
  rows.use <- rep(TRUE,nrow(models))
  cols.use <- Matrix::colSums(drop) > 0
  ## n <- unique(Matrix::rowSums(drop))
  n <- Matrix::rowSums(drop)
  ## index <- split(1:length(n),n)
  ## if(length(n)!=1)
  ##   stop("drop must be a model matrix with a fixed number of included snps (columns==1) per model (row).")
  drop <- drop[,cols.use]
  models.check <- models[,colnames(drop)]
  result <- drop %*% Matrix::t(models.check)
  ## system.time(to.drop <- lapply(index, function(i) {
  ##   Matrix::colSums(result[i,,drop=FALSE] == unique(n[i])) > 0
  ## }))
  to.drop <- colSums(result==n) > 0
  if(!quiet)
    cat(sum(to.drop),"identified (",sum(to.drop)/length(to.drop),").\n")
  models <- models[!to.drop,,drop=FALSE]
  if(!quiet)
    cat(nrow(models),"remain.\n")
  return(models)
}

##' Collapse a model matrix by SNP groups
##'
##' Internal function
##' @title models.group.collapse
##' @param models 
##' @param groups 
##' @return A dgCMatrix, with ncol==length(groups), and entry 1 if a SNP from that group is in the model.  Returns the input matrix in the case groups has zero length.
##' @author Chris Wallace
models.group.collapse <- function(models, groups) {
  if(length(groups)==0)
    return(models)
  newmodels <-  Matrix(0, nrow(models), length(groups),
                       dimnames=list(NULL,names(groups)),sparse=TRUE)
  for(i in 1:length(groups)) {
    if(length(groups[[i]])>1) {
      newmodels[,i] <- apply(models[, groups[[i]] ], 1, max)
    } else {
      newmodels[,i] <- models[, groups[[i]] ]
    }
  }
  return(newmodels)
}

models.add.excluded <- function(models, snps.exclude) {
  cbind2(models,
         Matrix(0, nrow(models), length(snps.exclude),
                dimnames=list(NULL,snps.exclude),sparse=TRUE))
}
 
snp.sort.models <- function(models, snps) {
  if(identical(snps, colnames(models)))
    return(models)
  
  return(models[,snps])
}

mexpand <- function(bma, groups) {
  models <- bma@models
  if(!all(names(groups) %in% colnames(models)))
    stop("All group index SNPs need to be in existing models\n")
  for(index.snp in names(groups)) {
    i <- which(models[,index.snp]==1)
    j <- which(colnames(models)==index.snp)
    msub <- models[i, -j]
    mkeep <- models[-i, -j]
    msub <- outer.models(models[i, -j, drop=FALSE],
                   make.models.single(groups[[index.snp]], n.use=1, quiet=TRUE))
    mkeep <- cbind2(models[-i, -j, drop=FALSE], Matrix(0, nrow(models)-length(i), length(groups[[index.snp]]),
                                          dimnames=list(NULL,groups[[index.snp]])))
    models <- rbind2(msub, mkeep)
  }
  return(models)
}

mgrow <- function(bma, quiet=FALSE) {

  models <- bma@models
  groups <- bma@groups
  snps <- colnames(models)
  snps.use <- snps[colSums(models)>0]
  snps.exclude <- setdiff(snps, snps.use)
  groups.use <- lapply(groups, setdiff, snps.exclude)
  models.new <- make.models(snps.use, n.use=1, groups=groups.use, quiet=TRUE)
  r.new <- nrow(models.new)
  r.parent <- nrow(models)
  rmodels.parent <- models[rep(1:r.parent, each=r.new),snps.use]
  rmodels.new <- models.new[rep(1:r.new, times=r.parent),snps.use]
  m <- rmodels.parent + rmodels.new

  ## ## prune1: SNPs should only be counted once
  ## system.time({rows.drop <- rowSums(m==2)>0
  ##              if(any(rows.drop))
  ##                m2 <- m[ !rows.drop, ]
  ##              mstring <- apply(m2,1,paste,collapse="")
  ##              m2 <- m2[ !duplicated(mstring), ]  
  ##  })

##  system.time({
   ## prune1: SNPs should only be counted once
  nzero <- summary(m)
  x2 <- which(nzero[,"x"]==2)
  i.drop <- nzero[ x2, "i"]
  nzero <- nzero[ !(nzero[,"i"] %in% i.drop), ]

  ## prune 2: multiple paths to each model, keep unique
  rows.split <- split( nzero[,"j"], nzero[,"i"] )
  rows.keep <- which(!duplicated(rows.split))
  m2 <- m[ as.numeric(names(rows.split)[rows.keep]), ]
  
  ## prune 3: max rs per group==1
  if(length(groups)) {
    rows.keep <- rep(TRUE,nrow(m2))
    for(g in groups) { ## UNTESTED!
      if(length(g)>1) {
        mg <- m2[,g]
        rs <- rowSums(mg)
        rows.keep <- rows.keep & rs <= 1
      }
    }
    m2 <- m2[ rows.keep, ]
  }

  if(length(snps.exclude))
    m2 <- models.add.excluded(m2, snps.exclude)

  return(snp.sort.models(m2, snps))

}

max.models.single <- function(n, n.use) {
  return(choose(n,n.use))
}

max.models.multi <- function(groups, n.use) {

  if(n.use==0)
    return(1)
  
  if(length(groups) < n.use)
    return(0)

  n.groups <- sapply(groups,length)
  
  combs <- combn(length(n.groups),n.use)
  n.models <- matrix( n.groups[ combs ], nrow=n.use )
  sum(apply( n.models, 2, prod ))
    
}

max.models <- function(snps, n.use, groups=list()) {

  if(!length(groups))
    return(max.models.single(length(snps),n.use))

  nogroups <- setdiff(unique(snps), unlist(groups))
  if(length(nogroups))
    groups <- c(groups, as.list(nogroups))
  n.groups <- sapply(groups, length)
  if(all(n.groups==1)) ## all singletons, ignore grouping
     return(max.models.single(length(snps),n.use))
  
  if(all(n.groups>1)) ## all multi, generate models using at most one SNP from each group
    return(max.models.multi(groups, n.use))

  ## mixed single and multi
  ## deal with single snp groups first - these are easy
  which.single <- which(n.groups==1)
  single.list <- sapply(0:n.use, function(i) max.models.single(length(which.single), i) )
  multi.list <- sapply(n.use:0, function(i) max.models.multi(groups[-which.single], i))  
  models.list <- single.list * multi.list
  return(sum(models.list))

}

## + should add colnames, rownames etc
