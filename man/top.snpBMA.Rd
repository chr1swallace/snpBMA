\name{top.snpBMA}
\alias{top.snpBMA}
\title{snpBMA}
\usage{
  top.snpBMA(object, nmodels = 6L, keep.snps = FALSE,
    what = "models")
}
\arguments{
  \item{object}{Object of class snpBMA}

  \item{nmodels}{number of models to show}

  \item{what}{\code{models} or \code{snps}, indicating
  whether to output the top models or SNPs}
}
\value{
  matrix of top models and their Bayes Factors
}
\description{
  Show top models in a snpBMA object
}
\author{
  Chris Wallace
}

