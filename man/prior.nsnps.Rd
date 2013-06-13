\name{prior.nsnps}
\alias{prior.nsnps}
\title{DELETE ME!}
\usage{
  prior.nsnps(x = 1:length(groups), groups, expected)
}
\arguments{
  \item{x}{number of SNPs in a model (defaults to
  1:length(groups), ie returns a vector)}

  \item{groups}{groups of SNPs, from which at most one SNP
  should be selected}

  \item{expected}{expected number of SNPs in a model}
}
\value{
  prior probability/ies as a numeric
}
\description{
  Binomial prior for number of SNPs in a model
}
\author{
  Chris Wallace
}

