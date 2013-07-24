\name{model.names}
\alias{model.names}
\title{model.names}
\usage{
  model.names(models, snps = NULL)
}
\arguments{
  \item{models}{model matrix}

  \item{num}{if TRUE (default), return numerics, otherwise
  return SNP names}
}
\value{
  character vector of length nrow(models)
}
\description{
  Make a character vector naming models uniquely according
  to the SNPs included
}
\author{
  Chris Wallace
}

