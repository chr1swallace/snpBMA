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
