##' Uses snpStats' snp.imputation function
##'
##' @title Fill in missing values in SnpMatrix object
##' @param X 
##' @param bp 
##' @param strata 
##' @param numeric 
##' @return matrix, either numeric (if numeric=TRUE) or SnpMatrix (if numeric=FALSE)
##' @author Chris Wallace
##' @export 
impute.missing <- function(X,bp=1:ncol(X),strata=NULL, numeric=FALSE) {
  N <- as(X,"numeric")
  if(!is.null(strata)) {
    strata <- as.factor(strata)
    if(length(levels(strata))>10)
      stop("too many levels in strata\n")
    for(i in levels(strata)) {
      cat("\nstrata",i,"\n")
      wh <- which(strata==i)
      N[wh,] <- impute.missing(X[wh,,drop=FALSE],bp, numeric=TRUE)
    }
  } else {
    csumm <- col.summary(X)
    wh <- which(csumm[,"Call.rate"]<1)
    cat(length(wh),"to impute\n")
    for(i in wh) {
      cat(i,".")
      rule <- snp.imputation(X[,-i,drop=FALSE],X[,i,drop=FALSE],bp[-i],bp[i])
      if(is.null(rule@.Data[[1]]))
        next
      imp <- impute.snps(rule,X[,-i,drop=FALSE])
      wh.na <- which(is.na(N[,i]))
      N[wh.na,i] <- imp[wh.na]
    }
    cat("\n")
  }
  if(numeric) {
    return(N)
  } else {
    return(new("SnpMatrix",
               data=round(N)+1, nrow=nrow(N), ncol=ncol(N), dimnames=dimnames(N)))
  }
}

