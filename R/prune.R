##' Compare potential parent and child models according to BF
##'
##' Returns models whose children SHOULD be visited
##' @title mcomp
##' @param bma snpBMA object
##' @param snps character vector of snps to drop
##' @return a snpBMA object with a subset of models containing only SNPs not included in \code{snps}
##' @author Chris Wallace
snps.prune.models <- function(bma, snps) {
  if(!all(snps %in% bma.1@snps))
    stop("SNPs to drop must be in BMA object\n")

  models <- bma@models
  models.drop <- which(rowSums(models[,snps]) > 0)
  if(length(models.drop))
    return(bma[-models.drop, ])

  return(bma)
}
 
##' Compare potential parent and child models according to BF
##'
##' Returns models whose children SHOULD be visited
##' @title mcomp
##' @param data snpBMAdata or snpBMAstrat object
##' @param snps character vector of snps to drop
##' @return object of same class as \code{data}, with SNPs in \code{snps} dropped.
##' @author Chris Wallace
snps.prune.data <- function(data, snps) {
  ## for each child model, identify its parents
  ## models to drop should be defined as the set with any 2*logbf(parent/child) > 2*log(rel.prior) + 2*lbf
  tags <- data@tags
  tags <- tags[!(tags %in% snps)]
  data <- data[,unique(tags)]
  data@tags <- tags
  return(data)    
##   make.data(data@.Data, data@Y,tags=tags, family=data@family, covar=data@covar,
##             strata=if(is(data,"snpBMAstrat")) { data@strata } else { NULL })

}
