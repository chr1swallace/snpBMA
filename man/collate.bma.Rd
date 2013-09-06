\name{collate.bma}
\alias{collate.bma}
\title{Collate a list of glib output}
\usage{
  collate.bma(results)
}
\arguments{
  \item{results}{a list of output from \code{\link{glib}}}
}
\value{
  an object similar to that returned by a single call of
  \code{glib}.  object$prior and object$posterior are _NOT_
  correct
}
\description{
  Internal: given a list of BMA results from glib, combine
  into a single object
}
\details{
  With the multicore package, it can be useful to farm out
  BMA jobs via mclapply
}
\author{
  Chris Wallace
}

