\name{null.model}
\alias{null.model}
\title{null.model}
\usage{
  null.model(snps = NULL, nsnps = length(snps))
}
\arguments{
  \item{snps}{character vector of snp names}

  \item{nsnps}{number of SNPs, defaults to length(snps)}
}
\value{
  Matrix, nrow=0, ncol=length(snps), colnames=snps
}
\description{
  Internal function, null model matrix
}
\details{
  Either snps or nsnps need to be supplied, not both.
}
\author{
  Chris Wallace
}

