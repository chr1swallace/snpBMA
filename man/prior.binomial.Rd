\name{prior.binomial}
\alias{prior.binomial}
\title{prior.binomial}
\usage{
  prior.binomial(x = 1:n, n = length(groups), expected,
    groups = NULL)
}
\arguments{
  \item{x}{number of SNPs in a model (defaults to
  1:length(groups), ie returns a vector)}

  \item{groups}{groups of SNPs, from which at most one SNP
  should be selected}

  \item{expected}{expected number of SNPs in a model}
}
\value{
  prior probability/ies as a numeric
}
\description{
  Binomial prior for number of SNPs in a model
}
\examples{
n<-100 # 100 SNPs in region
x <- 1:10 # consider prior for up to 10 causal SNPs
xbar <- 3 # expect around 3 causal

## a binomial prior
y <- prior.binomial(x, n, xbar)
plot(x, y, type="h")
}
\author{
  Chris Wallace
}

