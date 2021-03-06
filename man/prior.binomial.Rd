\name{prior.binomial}
\alias{prior.binomial}
\title{prior.binomial}
\usage{
  prior.binomial(x = 1:n, n = length(groups), expected,
    groups = NULL, pi0 = 0.2, value = c("odds", "prob"))
}
\arguments{
  \item{x}{number of SNPs in a model (defaults to
  1:length(groups), ie returns a vector)}

  \item{n}{total number of SNPs or SNP groups available}

  \item{expected}{expected number of SNPs in a model}

  \item{groups}{groups of SNPs, from which at most one SNP
  should be selected}

  \item{pi0}{prior probability that no SNP is associated}

  \item{value}{"odds" or "prob" for prior odds (relative to
  pi0) or prior probabilities}
}
\value{
  prior probability/ies or odds as a numeric
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

