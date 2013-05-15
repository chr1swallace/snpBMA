setClass("snpBMAdata",
         representation(X="matrix",
                        Y="numeric",
                        family="character",
                        tags="character"),
         validity=function(object) {
           if(length(object@Y) != nrow(object@X))
             stop("phenotype vector Y must have length == nrow(X)")
         })
##' 
setClass("snpBMA",
         representation(nsnps="numeric",
                        snps="character",
                        bf="matrix",
                        models="matrix",
                        glim.est="list"))
