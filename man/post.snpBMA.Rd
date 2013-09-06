\name{post.snpBMA}
\alias{post.snpBMA}
\title{Internal, post.snpBMA}
\usage{
  post.snpBMA(object, priors, snpsep = "/")
}
\arguments{
  \item{object}{Object of class snpBMA}

  \item{priors}{prior odds for models containing a fixed
  number of SNPs}

  \item{snpsep}{character string used to separate SNP names
  when naming models}
}
\value{
  pp matrix
}
\description{
  post.snpBMA
}
\details{
  Calculate posterior probs given BF and prior odds
}
\author{
  Chris Wallace
}

