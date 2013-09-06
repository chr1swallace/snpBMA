\name{top.snpBMAlist}
\alias{top.snpBMAlist}
\title{top.snpBMAlist}
\usage{
  top.snpBMAlist(object, priors, nmodels = 6L,
    min.pp = NULL, min.cpp = NULL, what = "models", ...)
}
\arguments{
  \item{object}{Object of class snpBMAlist}

  \item{priors}{vector of priors from which prior odds for
  models containing a fixed number of SNPs may be derived}

  \item{nmodels}{number of models to show}

  \item{min.pp}{("minimim posterior probability") optional,
  if present, limit to models with psterior probabilities >
  min.pp}

  \item{min.cpp}{("minimum cumulative posterior
  probability") optional, if present, limit to the best
  models which together account for at least min.cpp of the
  posterior probability}

  \item{what}{\code{models} or \code{snps}, indicating
  whether to output the top models or SNPs}

  \item{...}{additional arguments passed to
  \link{post.snpBMA}}
}
\value{
  matrix of top models and their Bayes Factors
}
\description{
  Show top models in a snpBMAlist object
}
\author{
  Chris Wallace
}

