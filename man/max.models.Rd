\name{max.models}
\alias{max.models}
\title{max.models}
\usage{
  max.models(snps = character(0), n.use,
    n.snps = length(snps), groups = list())
}
\arguments{
  \item{snps}{Character vector of snp names covering all
  snps in the region.  Needed if you want to supply
  \code{groups}}

  \item{n.use}{Number of SNPs in model}

  \item{n.snps}{Number of SNPs in region, length(snps).
  Either \code{snps} or \code{n.snps} must be specified.}

  \item{groups}{optional list of character vectors, each
  giving groups of SNPs in LD, so that only one of any
  group is included in a model}
}
\value{
  The number of possible models
}
\description{
  Total possible models that may be formed of n.use out of
  snps SNPs
}
\details{
  This is a choose(snps, n.use) in the simple case, but
  subject to the groups list.Total possible models
}
\author{
  Chris Wallace
}

