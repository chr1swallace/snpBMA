\name{snp.summary}
\alias{snp.summary}
\title{Summary Bayes Factors for each SNP}
\usage{
  snp.summary(object, do.order = TRUE)
}
\arguments{
  \item{object}{Object of class snpBMA}

  \item{do.order}{if TRUE (the default), reorder the
  returned summary of snp Bayes Factors, so that those with
  the strongest support are shown first}
}
\value{
  a Matrix of class \code{\linkS4class{dgeMatrix}} with
  each row corresponding to a SNP, and each column to a
  different value of phi
}
\description{
  Summary Bayes Factors for each SNP
}
\details{
  The summary Bayes Factor is the mean of the Bayes Factor
  for all models containing that SNP, which is valid if
  each model has equal prior
}
\author{
  Chris Wallace
}

