################################################################################
## stack
##' @rdname stack-methods
##' @aliases stack,snpBMA-method
setMethod("stack", signature(object="snpBMA"),
          function(object, ...) {
              x <- list(object, ...)
              null <- sapply(x,is.null)
              new("snpBMAlist", x[!null]) })
##' @rdname stack-methods
##' @aliases stack,modelSummary-method 
setMethod("stack",signature(object="modelSummary"),
          function(object, ...) {
            L <- list(object,...)
            x <- do.call("rbind", L)
            colnames(x) <- colnames(object)
            snps <- do.call("c",lapply(L, function(o) o@snps))
            snpsep <- object@snpsep
            new("modelSummary",
                .Data=x, snps=snps, snpsep=snpsep)
          })

################################################################################
## show

##' Show snpBMA and related objects
##'
##' For all objects, a very brief description of the object size is
##' given.  For snpBMA objects only, the models with greatest support are shown.
##' @title show-methods
##' @param object object of class snpBMA, snpBMAdata or snpBMAlist
##' @return no return value
##' @author Chris Wallace
##' @export
##' @docType methods
##' @keywords methods
##' @rdname show-methods
##' @aliases show,snpBMA-method
setMethod("show", signature="snpBMA",
          function(object) {
            nmod <- nrow(object@models)
            nsnp <- ncol(object@models)
            n.use <- apply(object@models,1,sum)
            cat("Bayes factors for",nmod,"models from",object@nmodels,"possible models using",nsnp,"SNPs.\nSNPs per model:\n")
            print(table(n.use))
            cat("Top models:\n")
            print(top.models(object))
          })
##' @rdname show-methods
##' @aliases show,snpBMAlist-method
setMethod("show", signature(object="snpBMAlist"),
          function(object) {
            cat("A list of",length(object),"snpBMA objects.\n") })
##' @rdname show-methods
##' @aliases show,modelSummary-method
setMethod("show", signature(object="modelSummary"),
          function(object) {
            tmp <- object@.Data
            rownames(tmp) <- lapply(object@snps, paste, collapse=object@snpsep)
            show(tmp) })
##' @rdname show-methods
##' @aliases show,snpBMAdata-method
setMethod("show", signature="snpBMAdata",
          function(object) {
            nsamp <- length(object@Y)
            nsnp <- length(object@tags)
            ntag <- ncol(object@.Data)
            cat("snpBMAdata object, with",object@family,"phenotypes on\n",nsamp,"individuals\n",nsnp,"SNPs, represented by\n",ntag,"tags.\n")
          })
##' @rdname show-methods
##' @aliases show,snpBMAstrat-method
setMethod("show", signature="snpBMAstrat",
          function(object) {
            nsamp <- length(object@Y)
            nsnp <- ncol(object@.Data)
            nstrat <- length(unique(object@strata))
            
            cat(class(object),"object, with",object@family,"phenotypes.\n",nsamp,"individuals in",nstrat,"strata;",nsnp,"SNPs represented by",length(unique(object@tags)),"tags.\n")
          })

################################################################################

setMethod("initialize", signature(.Object = "snpBMA"),
    function (.Object, nmodels, snps, groups, bf, models, nsnps, ...)  {
      if(length(groups)) {
        nogroup <- setdiff(unique(snps), unlist(groups))
        names(nogroup) <- nogroup
        if(length(nogroup))          
          groups <- c(groups, as.list(nogroup))
      }
      callNextMethod(.Object, nmodels=nmodels, snps=snps, groups=groups, bf=bf,
                     models=models, nsnps=nsnps, ...)
    })
setMethod("[",
          signature=c(x="snpBMA", i="ANY", j="missing", drop="missing"),
          function(x, i) {
            new("snpBMA",
              nmodels=x@nmodels,
              snps=x@snps,
              groups=x@groups,
              bf=x@bf[i, , drop=FALSE ],
              models=x@models[i, , drop=FALSE ],
              nsnps=x@nsnps) })
setMethod("[",
          signature=c(x="modelSummary", i="ANY", j="missing", drop="ANY"),
          function(x, i) {
            new("modelSummary",
                .Data=x@.Data[i,,drop=FALSE],
                snps=x@snps[i],
                snpsep=x@snpsep)})

setMethod("snps0",
          signature=c(object="snpBMA"),
          function(object) {
            snps <- colnames(object@models)
            names(snps) <- 0:(length(snps)-1)
            return(snps)
          })

################################################################################
## snps

#' @rdname snps-methods
#' @aliases snps,modelSummary-method
setMethod("snps",
          signature=c(object="modelSummary"),
          function(object) {
            unique(unlist(object@snps))                   
          })

################################################################################
## top.models

#' @rdname top.models-methods
#' @aliases top.models,snpBMA-method
#' @seealso \code{\link{top.snpBMA}} for application to top.snpBMA objects
setMethod("top.models",
          signature=c(object="snpBMA"),
          function(object, ...) {
            top.snpBMA(object, ...) })

#' @rdname top.models-methods
#' @aliases top.models,snpBMAlist-method
#' @seealso \code{\link{top.snpBMAlist}} for application to top.snpBMAlist objects
setMethod("top.models",
          signature=c(object="snpBMAlist"),
          function(object, ...) {
            top.snpBMAlist(object, ...) })

################################################################################
## top.snps

#' @rdname top.snps-methods
#' @aliases top.snps,snpBMA-method
setMethod("top.snps",
          signature=c(object="snpBMA"),
          function(object, ...) {
            top.snpBMA(object, what="snps", ...) })

#' @rdname top.snps-methods
#' @aliases top.snps,snpBMAlist-method
setMethod("top.snps",
          signature=c(object="snpBMAlist"),
          function(object, ...) {
            top.snpBMAlist(object, what="snps", ...) })



################################################################################

## snpBMAlist
## setMethod("initialize", signature(.Object = "snpBMAlist"),
##     function (.Object, ...)  {
##       .Object@object=list(...)
##       return(.Object)  })
################################################################################
##  snpBMAdata
setMethod("[",
          signature=c(x="snpBMAdata", i="missing", j="ANY", drop="missing"),
          function(x, j) {
            newX <- x@.Data[,j,drop=FALSE]
            newsnps <- colnames(newX)
            newtags <- x@tags[ x@tags %in% newsnps ]
            new("snpBMAdata",
                .Data=newX,
                Y=x@Y,
                family=x@family,
                tags=newtags,
                covar=x@covar) })
setMethod("[",
          signature=c(x="snpBMAdata", i="ANY", j="missing", drop="missing"),
          function(x, i) {
            covar <- x@covar
            if(nrow(covar))
              covar <- covar[i,,drop=FALSE]
            new("snpBMAdata",
                .Data=x@.Data[i,,drop=FALSE],
                Y=x@Y[i],
                family=x@family,
                tags=x@tags,
                covar=covar) })
 setMethod("[",
          signature=c(x="snpBMAstrat", i="missing", j="ANY", drop="missing"),
          function(x, j) {
            newX <- x@.Data[,j,drop=FALSE]
            newsnps <- colnames(newX)
            newtags <- x@tags[ x@tags %in% newsnps ]            
            new("snpBMAstrat",
                .Data=newX,
                Y=x@Y,
                family=x@family,
                tags=newtags,
                covar=x@covar,
                strata=x@strata) })
setMethod("[",
          signature=c(x="snpBMAstrat", i="ANY", j="missing", drop="missing"),
          function(x, i) {
            newstrat=x@strata[i]
            covar <- x@covar
            if(nrow(covar))
              covar <- covar[i,,drop=FALSE]
          if(length(unique(newstrat))==1) {
              new("snpBMAdata",
                  .Data=x@.Data[i,,drop=FALSE],
                  Y=x@Y[i],
                  family=x@family,
                  tags=x@tags,
                  covar=covar)
            } else {
              new("snpBMAstrat",
                  .Data=x@.Data[i,,drop=FALSE],
                  Y=x@Y[i],
                  family=x@family,
                  tags=x@tags,
                  covar=covar,
                  strata=newstrat)
            } })

################################################################################
## Models
setMethod("show", signature="Models",
          function(object) {
            nmod <- nrow(object@models)
            nsnp <- ncol(object@models)
            cat(nmod,"models, each including",object@nsnps,"out of",nsnp,"SNPs.\n")
          })
setMethod("models","Models",
          function(object) {
            return(object@models)
          })

setClass("dropModels",
         representation(snps="character"),
         contains="Models",
         validity=function(object) {
           if(object@snps %in% colnames(object@models))
             stop("models to be dropped must be formed only from SNPs which have NOT been dropped.")
         })

################################################################################
## dropModels
setMethod("initialize", signature(.Object = "dropModels"),
    function (.Object, models, snps, ...)  {
      if(any(snps %in% colnames(models))) {
        ntarget <- unique(Matrix::rowSums(models))
        if(length(ntarget)>1)
          stop("models matrix should have a fixed number of snps in each model, but varies")
        models <- models[,setdiff(colnames(models),snps),drop=FALSE]
        if(length(ntarget)==1) {
          rs <- Matrix::rowSums(models)
          models <- models[rs==ntarget,,drop=FALSE]
        }        
      }
      .Object@models <- models
      .Object@snps <- snps
      return(.Object)
##      callNextMethod(.Object, models=models, snps=snps, ...)
    })
setMethod("show", signature(object = "dropModels"),
    function (object) {
      cat("Models to drop:\t", nrow(object@models),
          "\nSNPs to drop:\t", length(object@snps),"\n") 
    })
setMethod("snps","dropModels",
          function(object) {
            return(object@snps)
          })
setMethod("stack","list",
          function(object) {
            nobj <- length(object)
            classes.ok <- sapply(object, is, "dropModels")
            if(!all(classes.ok))
              stop("object must be a list of objects of class dropModels.")
            snps.drop <- unique(unlist(lapply(object, snps)))
            object <- lapply(object, function(obj) {
              new("dropModels", models=obj@models, snps=snps.drop)
            })
            new("dropModels",
                models=do.call("rBind",lapply(object, models)),
                snps=snps.drop)
          })

################################################################################
## prune
setMethod("prune",
    signature(object = "snpBMA", prune.object = "snpBMA"),
    function (object, prune.object, ...) 
    {
      models.prune(children=object, parents=prune.object, ...)
    }
)
##' @rdname snps.prune-methods
##' @aliases snps.prune,snpBMA-method
setMethod("snps.prune",
    signature(object = "snpBMA", snps = "character"),
    function (object, snps, ...) 
    {
      snps.prune.models(bma=object, snps=snps)
    }
)
##' @rdname snps.prune-methods
##' @aliases snps.prune,snpBMAdata-method
setMethod("snps.prune",
    signature(object = "snpBMAdata", snps = "character"),
    function (object, snps, ...) 
    {
      snps.prune.data(data=object, snps=snps)
    }
)
##' @rdname snps.prune-methods
##' @aliases snps.prune,snpBMAstrat-method
setMethod("snps.prune",
    signature(object = "snpBMAstrat", snps = "character"),
    function (object, snps, ...) 
    {
      snps.prune.data(data=object, snps=snps)
    }
)
