\name{null.model}
\alias{null.model}
\title{null.model}
\usage{
  null.model(snps = NULL, nsnps = length(snps))
}
\arguments{
  \item{snps}{character vector of snp names}
}
\value{
  Matrix, nrow=0, ncol=length(snps), colnames=snps
}
\description{
  Internal function, null model matrix
}
\author{
  Chris Wallace
}

