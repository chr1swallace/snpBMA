\name{models.prune.snps}
\alias{models.prune.snps}
\title{mcomp}
\usage{
  models.prune.snps(bma, snps.to.drop, ...)
}
\arguments{
  \item{bma}{snpBMA object}

  \item{snps}{character vector of snps to drop}

  \item{...}{arguments passed to mcomp.detail()}
}
\value{
  a snpBMA object with a subset of models in the input
  \code{children} object, all child models for which
  twologBF (parent/child) > 2*log(pp.fold)
}
\description{
  Compare potential parent and child models according to BF
}
\details{
  Returns models whose children SHOULD be visited
}
\author{
  Chris Wallace
}

