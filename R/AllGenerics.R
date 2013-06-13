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

##' Stack snpBMA objects together
##'
##' @title stack-methods
##' @param object object of class snpBMA
##' @param ... more objects of class snpBMA
##' @return object of class snpBMAlist
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
##' @param ... number of models to show
##' @return no return value, models are printed to screen
##' @author Chris Wallace
##' @export
##' @docType methods
##' @rdname top.models-methods
##' @aliases top.models-methods
##' @keywords methods
setGeneric("top.models", function(object, ...) { standardGeneric("top.models") })

setGeneric("prune", function(object, prune.object, ...) { standardGeneric( "prune" ) })

