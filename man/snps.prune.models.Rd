\name{snps.prune.models}
\alias{snps.prune.models}
\title{mcomp}
\usage{
  snps.prune.models(bma, snps)
}
\arguments{
  \item{bma}{snpBMA object}

  \item{snps}{character vector of snps to drop}
}
\value{
  a snpBMA object with a subset of models containing only
  SNPs not included in \code{snps}
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

