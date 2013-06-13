\name{bma.nsnps}
\alias{bma.nsnps}
\title{bma.nsnps}
\usage{
  bma.nsnps(data, nsnps = 1, groups = list(),
    models.drop = NULL)
}
\arguments{
  \item{data}{object of class snpBMAdata}

  \item{nsnps}{number of SNPs to include in each model}

  \item{groups}{list of character vectors.  Each vector
  contains a set of SNPs in LD, of which at most one should
  be included in any model.}

  \item{models.drop}{matrix of models that shouldn't be
  visited.  This option is deprecated, it is better use
  \code{\link{bma.grow}()} instead.}
}
\value{
  object of class snpBMA
}
\description{
  Fit possible models with a fixed number of SNPs in each
  model
}
\author{
  Chris Wallace
}

