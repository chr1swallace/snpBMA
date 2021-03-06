\name{logdiff}
\alias{logdiff}
\title{logdiff (copied from coloc package)}
\usage{
  logdiff(x, y)
}
\arguments{
  \item{x}{numeric}

  \item{y}{numeric}
}
\value{
  max(x) + log(exp(x - max(x,y)) - exp(y-max(x,y)))
}
\description{
  Internal function, logdiff
}
\details{
  This function calculates the log of the difference of the
  exponentiated logs taking out the max, i.e. insuring that
  the difference is not negative
}
\author{
  Chris Wallace
}

