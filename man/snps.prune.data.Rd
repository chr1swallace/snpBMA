\name{snps.prune.data}
\alias{snps.prune.data}
\title{mcomp}
\usage{
  snps.prune.data(data, snps)
}
\arguments{
  \item{data}{snpBMAdata or snpBMAstrat object}

  \item{snps}{character vector of snps to drop}
}
\value{
  object of same class as \code{data}, with SNPs in
  \code{snps} dropped.
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

