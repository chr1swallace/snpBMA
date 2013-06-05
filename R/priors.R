##' Binomial prior for number of SNPs in a model
##'
##' @title prior.nsnps
##' @param x number of SNPs in a model (defaults to 1:length(groups), ie returns a vector)
##' @param groups groups of SNPs, from which at most one SNP should be selected
##' @param expected expected number of SNPs in a model
##' @return prior probability/ies as a numeric
##' @export
##' @author Chris Wallace
##' @examples
##' n<-100 # 100 SNPs in region
##' x <- 1:10 # consider prior for up to 10 causal SNPs
##' xbar <- 3 # expect around 3 causal
##'
##' ## a binomial prior
##' y <- prior.binomial(x, n, xbar)
##' plot(x, y, type="h")
prior.binomial <- function(x=1:n, n=length(groups), expected, groups=NULL) {
  x <- as.integer(x)
  if(any(x<0))
    stop("x should be an integer vector >= 0")
  if(any(x>n))
    stop("max x should be <= n")
  p <- expected/n
  dbinom(x, size=n, prob=p)  
}

##' Beta Binomial prior for number of SNPs in a model
##'
##' A binomial prior for the number of SNPs in a model may be
##' considered too peaked if there is relatively little prior
##' information about the number of causal SNPs, and, particularly if
##' the posterior model choice is sensitive to this prior, it can be
##' useful to consider a prior with greater spread.  One such choice
##' is the beta binomial model, implemented here.
##' 
##' @title prior.nsnps
##' @param x number of SNPs in a model (defaults to 1:length(groups), ie returns a vector)
##' @param groups groups of SNPs, from which at most one SNP should be selected
##' @param expected expected number of SNPs in a model
##' @param overdispersion overdispersion parameter.  Setting this to 1
##' gives a binomial prior.  Values < 1 are nonsensical: if you really
##' believe the prior should be underdispersed relative to a binomial
##' distribution, consider using a hypergeometric prior.
##' @return prior probability/ies as a numeric vector
##' @export
##' @examples
##' n<-100 # 100 SNPs in region
##' x <- 1:10 # consider prior for up to 10 causal SNPs
##' xbar <- 3 # expect around 3 causal
##'
##' ## a binomial prior
##' y <- prior.binomial(x, n, xbar)
##' plot(x, y, type="h")
##' 
##' ## is equivalent to
##' y1.0 <- prior.betabinomial(x, n, xbar, overdispersion=1.0)
##' points(x, y1.0, col="red")
##'
##' ##larger values of overdispersion change the distribution:
##' y1.1 <- prior.betabinomial(x, n, xbar, overdispersion=1.1)
##' y1.5 <- prior.betabinomial(x, n, xbar, overdispersion=1.5)
##' y2.0 <- prior.betabinomial(x, n, xbar, overdispersion=2.0)
##' points(x, y1.1, col="orange")
##' points(x, y1.5, col="pink")
##' points(x, y2.0, col="green")
##' @author Chris Wallace
prior.betabinomial <- function(x=1:n, n=length(groups), expected, overdispersion=1, groups=NULL) {
  if(overdispersion < 1)
    stop("overdispersion parameter should be >= 1")
   x <- as.integer(x)
  if(any(x<0))
    stop("x should be an integer vector >= 0")
  if(any(x>n))
    stop("max x should be <= n")
  p <- expected/n
  rho <- (overdispersion - 1)/(n-1)
  dbetabinom(x, size=n, prob=p, rho=rho)             
}

solve.bb <- function(pi, rho) {
  alpha <- (1-rho)*pi/(rho*pi + rho*(1-pi))
  beta <- alpha*(1-pi)/pi
  return(list(alpha=alpha,beta=beta))
}
