\name{make.data}
\alias{make.data}
\title{make.data}
\usage{
  make.data(X, Y, tags, family = "binomial")
}
\arguments{
  \item{X}{matrix of explanatory variables, either a
  \code{SnpMatrix} object, or a numeric matrix with 0, 1, 2
  indicating copies of a counted allele at each SNP.  Rows
  are samples, columns are SNPs.}

  \item{Y}{response}

  \item{tags}{a named character vector created by the
  \code{\link{tag}()} function}

  \item{family}{character string, "binomial" or "gaussian",
  passed to glib}
}
\value{
  a list with entries X and Y
}
\description{
  Create an snpBMAdata object given a matrix of genotypes
  and a vector of phenotypes
}
\author{
  Chris Wallace
}

