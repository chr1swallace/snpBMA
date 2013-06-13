\name{prior.betabinomial}
\alias{prior.betabinomial}
\title{prior.betabinomial}
\usage{
  prior.betabinomial(x = 1:n, n = length(groups), expected,
    overdispersion = 1, groups = NULL)
}
\arguments{
  \item{x}{number of SNPs in a model (defaults to
  1:length(groups), ie returns a vector)}

  \item{groups}{groups of SNPs, from which at most one SNP
  should be selected}

  \item{expected}{expected number of SNPs in a model}

  \item{overdispersion}{overdispersion parameter.  Setting
  this to 1 gives a binomial prior.  Values < 1 are
  nonsensical: if you really believe the prior should be
  underdispersed relative to a binomial distribution,
  consider using a hypergeometric prior.}
}
\value{
  prior probability/ies as a numeric vector
}
\description{
  Beta Binomial prior for number of SNPs in a model
}
\details{
  A binomial prior for the number of SNPs in a model may be
  considered too peaked if there is relatively little prior
  information about the number of causal SNPs, and,
  particularly if the posterior model choice is sensitive
  to this prior, it can be useful to consider a prior with
  greater spread.  One such choice is the beta binomial
  model, implemented here, under which the number of SNPs
  follows a binomial distribution with parameters n, p
  while p follows a beta distribution with parameters
  chosen so that the mean and the overdispersion (relative
  to a binomial distribution) of the number of SNPs is as
  specified by \code{expected} and \code{overdispersion},
  respectively.
}
\examples{
n<-100 # 100 SNPs in region
x <- 1:10 # consider prior for up to 10 causal SNPs
xbar <- 3 # expect around 3 causal

## a binomial prior
y <- prior.binomial(x, n, xbar)
plot(x, y, type="h")

## is equivalent to
y1.0 <- prior.betabinomial(x, n, xbar, overdispersion=1.0)
points(x, y1.0, col="red")

##larger values of overdispersion change the distribution:
y1.1 <- prior.betabinomial(x, n, xbar, overdispersion=1.1)
y1.5 <- prior.betabinomial(x, n, xbar, overdispersion=1.5)
y2.0 <- prior.betabinomial(x, n, xbar, overdispersion=2.0)
points(x, y1.1, col="orange")
points(x, y1.5, col="pink")
points(x, y2.0, col="green")
}
\author{
  Chris Wallace
}

