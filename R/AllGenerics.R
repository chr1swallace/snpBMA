##' Extract model matrix
##'
##' @title models
##' @param object object of class snpBMA 
##' @return matrix
##' @author Chris Wallace
##' @export
##' @docType methods
##' @rdname models-methods
setGeneric("models", function(object){ standardGeneric ("models") })

setGeneric("snps", function(object){ standardGeneric ("snps") })

setGeneric("snps0", function(object){ standardGeneric ("snps0") })

setGeneric("stackModels", function(object) { standardGeneric ("stackModels") })

##' Stack objects together
##'
##' @title stack-methods
##' @param object object of class snpBMA or modelSummary
##' @param ... more objects of the same class
##' @return object of class snpBMAlist or modelSummary
##' @author Chris Wallace
##' @export
##' @docType methods
##' @rdname stack-methods
##' @aliases stack-methods
##' @keywords methods
setGeneric("stack", function(object,...) { standardGeneric ("stack") })

##' Show models with greatest posterior probability
##'
##' @title top.models
##' @param object object of class snpBMA or snpBMAlist
##' @param ... additional arguments, notably \code{nmodels}: number of models to show
##' @return no return value, models are printed to screen
##' @author Chris Wallace
##' @export
##' @docType methods
##' @rdname top.models-methods
##' @aliases top.models-methods
##' @keywords methods
setGeneric("top.models", function(object, ...) { standardGeneric("top.models") })

setGeneric("prune", function(object, prune.object, ...) { standardGeneric( "prune" ) })

##' Drop SNPs from a snpBMAdata or snpBMAstrat object
##'
##' @title snps.prune
##' @param object a snpBMAdata or snpBMAstrat object
##' @param prune.snps character vector of SNPs to drop
##' @param ... other arguments, currently ignored
##' @return object of the same class as \code{object}
##' @export
##' @author Chris Wallace
##' @docType methods
##' @rdname snps.prune-methods
##' @aliases snps.prune-methods
##' @keywords methods
setGeneric("snps.prune", function(object, snps, ...) { standardGeneric( "snps.prune" ) })
