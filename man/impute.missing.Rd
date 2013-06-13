\name{impute.missing}
\alias{impute.missing}
\title{Fill in missing values in SnpMatrix object}
\usage{
  impute.missing(X, bp = 1:ncol(X), strata = NULL,
    numeric = FALSE)
}
\arguments{
  \item{X}{}

  \item{bp}{}

  \item{strata}{}

  \item{numeric}{}
}
\value{
  matrix, either numeric (if numeric=TRUE) or SnpMatrix (if
  numeric=FALSE)
}
\description{
  Uses snpStats' snp.imputation function
}
\author{
  Chris Wallace
}

