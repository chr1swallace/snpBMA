\name{make.models.single}
\alias{make.models.single}
\title{make.models.single}
\usage{
  make.models.single(snps, n.use, quiet = FALSE)
}
\arguments{
  \item{snps}{character vector of snp names}

  \item{n.use}{integer: number of snps to use in each
  model}

  \item{quiet}{if TRUE, suppress progress messages}
}
\value{
  Matrix, with nrow=number of models, ncol=number of SNPs,
  entries of 1/0 indicate whether SNP is in/excluded from
  model
}
\description{
  Internal function, return a model matrix for choosing
  n.use snps out of snps
}
\author{
  Chris Wallace
}

