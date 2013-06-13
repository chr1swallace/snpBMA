\name{show.ld}
\alias{show.ld}
\title{show.ld}
\usage{
  show.ld(snps = colnames(X), samples = rownames(X),
    X = snp.data, lines.limit = 20)
}
\arguments{
  \item{snps}{optional character vector of column names of
  X (SNPs) for which LD should be calculated}

  \item{samples}{optional character vector of row names of
  X (Samples) for which LD should be calculated}

  \item{X}{a SnpMatrix object}

  \item{lines.limit}{}
}
\value{
  invisibly returns matrix with D' in upper.tri() entries
  and R^2 in lower.tri() entries
}
\description{
  Calculate D' and R^2 for a SnpMatrix object and display
  as a heatmap
}
\details{
  Side effect: displays plot on current graphics device
}
\author{
  Chris Wallace
}

