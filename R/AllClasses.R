################################################################################
## Models

setClass("Models",
         representation(models="dgCMatrix",
                        nsnps="numeric"),
         validity=function(object) {
           if(!all(rowSums(object@models) == object@nsnps))
             stop("Each model matrix should contain models with a fixed colSums()")
         })

## Model summary matrix
setClass("modelSummary",
         representation(snps="list",
                        snpsep="character"),
         contains="matrix",
         validity=function(object) {
           if(nrow(object@.Data)!=length(object@snps))             
             stop("snps list must index rows of .Data matrix")
           if(length(object@snpsep)!=1)
             stop("snpsep must be a single character")
         })

################################################################################
## snpBMAdata
setClass("snpBMAdata",
         representation(Y="numeric",
                        family="character",
                        tags="character",
                        covar="matrix"),
         contains="matrix",
         validity=function(object) {
           if(length(object@Y) != nrow(object@.Data))
             stop("phenotype vector Y must have length == nrow(X)")
           if(nrow(object@covar) > 0 && nrow(object@covar) != length(object@Y))
             stop("covar variables must have length == nrow(X)")
         })

setClass("snpBMAstrat",
         representation(.Data="list"),
         validity=function(object) {
           classes <- lapply(object@.Data,class)
           if(!all(classes=="snpBMAdata"))
             stop("all elements in a snpBMAstrat object must be of class snpBMAdata")
         })

################################################################################
## snpBMA
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
