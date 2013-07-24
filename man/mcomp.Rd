\name{mcomp}
\alias{mcomp}
\title{mcomp}
\usage{
  mcomp(parents, children, ...)
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
  .. content for \details{} ..
}
\author{
  Chris Wallace
}

