\name{models.prune}
\alias{models.prune}
\title{mcomp}
\usage{
  models.prune(parents, children, ...)
}
\arguments{
  \item{parents}{snpBMA object}

  \item{children}{snpBMA object}

  \item{...}{arguments passed to mcomp.detail()}
}
\value{
  a snpBMA object with a subset of models in the input
  \code{children} object, all child models for which
  twologBF (parent/child) > 2*log(pp.fold)
}
\description{
  Compare potential parent and child models according to BF
}
\details{
  Returns models whose children SHOULD be visited
}
\author{
  Chris Wallace
}

