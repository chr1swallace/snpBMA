\name{make.models.multi}
\alias{make.models.multi}
\title{make.models.multi}
\usage{
  make.models.multi(snps, n.use, groups, quiet = FALSE)
}
\arguments{
  \item{snps}{character vector of snp names}

  \item{n.use}{integer: number of snps to use in each
  model}

  \item{groups}{list of character vectors, each vector
  indicating a group of SNPs from which a maximum of one
  may be selected}

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

