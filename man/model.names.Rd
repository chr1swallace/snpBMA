\name{model.names}
\alias{model.names}
\title{model.names}
\usage{
  model.names(models, snps = NULL)
}
\arguments{
  \item{models}{model matrix}

  \item{snps}{optional, if present, restrict by snps}
}
\value{
  character vector of length nrow(models)
}
\description{
  Internal: make a character vector naming models uniquely
  according to the SNPs included
}
\author{
  Chris Wallace
}

