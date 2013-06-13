\name{models.diff}
\alias{models.diff}
\title{mcomp}
\usage{
  models.diff(parents, children, ...)
}
\arguments{
  \item{parents}{snpBMA object}

  \item{children}{snpBMA object}

  \item{...}{arguments passed to mcomp.detail()}
}
\value{
  a Matrix containing all child models for which twologBF
  (parent/child) > 2*log(rel.prior)
}
\description{
  Compare potential parent and child models according to BF
}
\details{
  Returns models whose children SHOULD NOT be visited
}
\author{
  Chris Wallace
}

