################################################################################
## Models

setClass("Models",
         representation(models="dgCMatrix",
                        nsnps="numeric"),
         validity=function(object) {
           if(!all(rowSums(object@models) == object@nsnps))
             stop("Each model matrix should contain models with a fixed colSums()")
         })
setMethod("show", signature="Models",
          function(object) {
            nmod <- nrow(object@models)
            nsnp <- ncol(object@models)
            cat(nmod,"models, each including",object@nsnps,"out of",nsnp,"SNPs.\n")
          })
setGeneric("models", function(object){ standardGeneric ("models") })
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
setGeneric("snps", function(object){ standardGeneric ("snps") })
setMethod("snps","dropModels",
          function(object) {
            return(object@snps)
          })
setGeneric("stackModels", function(object) { standardGeneric ("stackModels") })
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
## snpBMA, snpBMAdata
setClass("snpBMAdata",
         representation(X="matrix",
                        Y="numeric",
                        family="character",
                        tags="character"),
         validity=function(object) {
           if(length(object@Y) != nrow(object@X))
             stop("phenotype vector Y must have length == nrow(X)")
         })
setClass("snpBMA",
         representation(nmodels="numeric",
                        snps="character",
                        groups="list",
                        bf="matrix"),
         contains="Models",
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
            cat("Bayes factors for",nmod,"models from",object@nmodels,"possible models using",nsnp,"SNPs.\nSNPs per model:\n")
            print(table(n.use))
            cat("Top models:\n")
            print(top.models(object))
          })

