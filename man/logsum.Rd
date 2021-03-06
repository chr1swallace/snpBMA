\name{logsum}
\alias{logsum}
\title{logsum (copied from coloc package)}
\usage{
  logsum(x)
}
\arguments{
  \item{x}{numeric vector}
}
\value{
  max(x) + log(sum(exp(x - max(x))))
}
\description{
  Internal function, logsum
}
\details{
  This function calculates the log of the sum of the
  exponentiated logs taking out the max, i.e. insuring that
  the sum is not Inf
}
\author{
  Claudia Giambartolomei
}

