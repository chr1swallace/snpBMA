################################################################################
## Models

setClass("Models",
         representation(models="dgCMatrix",
                        nsnps="numeric"),
         validity=function(object) {
           if(!all(rowSums(object@models) == object@nsnps))
             stop("Each model matrix should contain models with a fixed colSums()")
         })

################################################################################
## snpBMA, snpBMAdata
setClass("snpBMAdata",
         representation(Y="numeric",
                        family="character",
                        tags="character"),
         contains="matrix",
         validity=function(object) {
           if(length(object@Y) != nrow(object@.Data))
             stop("phenotype vector Y must have length == nrow(X)")
         })

setClass("snpBMAstrat",
         representation(strata="factor"),
         contains="snpBMAdata",
         validity=function(object) {
           if(length(object@Y) != nrow(object@.Data))
             stop("phenotype vector Y must have length == nrow(X)")
          if(length(object@strata) != nrow(object@.Data))
             stop("strata vector must have length == nrow(X)")
         })

setClass("snpBMA",
         representation(nmodels="numeric",
                        snps="character",
                        groups="list",
                        bf="matrix",
                        prior="list"),
         contains="Models",
         validity=function(object) {
           if(nrow(object@bf) != nrow(object@models))
             stop("bf must have same nrow as models")
           if(length(unique(object@snps)) != ncol(object@models))
             stop("models must have ncol == length(unique(snps))")
         })

setClass("snpBMAlist",
         contains="list",
         validity=function(object) {
           if(!all(sapply(object, is, "snpBMA")))
             stop("object must be a list of snpBMA objects.\n")
         })
