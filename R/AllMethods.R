################################################################################
## snpBMA, snpBMAdata
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

setMethod("[",
          signature=c(x="snpBMA", i="ANY", j="missing", drop="missing"),
          function(x, i) {
            new("snpBMA",
              nmodels=x@nmodels,
              snps=x@snps,
              groups=x@groups,
              bf=x@bf[i, ],
              models=x@models[i, ],
              nsnps=x@nsnps) })
  

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
setMethod("stackModels","list",
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
setMethod("prune",
    signature(object = "snpBMA", prune.object = "character"),
    function (object, prune.object) 
    {
      models.prune.snps(children=object, parents=prune.object, ...)
    }
)
