\name{mcomp.detail}
\alias{mcomp.detail}
\title{mcomp.detail}
\usage{
  mcomp.detail(m.parent, m.child, bf.parent, bf.child,
    ntarget, n.child, prior.parents, prior.children, what,
    pp.fold = 10, quiet = FALSE)
}
\arguments{
  \item{m.parent}{parent model matrix}

  \item{m.child}{child model matrix}

  \item{bf.parent}{BF matrix for parent models}

  \item{bf.child}{BF matrix for child models}

  \item{ntarget}{number of SNPs in parent models}

  \item{n.child}{number of SNPs in child models}

  \item{prior.parents}{prior for parent models}

  \item{prior.children}{prior for child models,}

  \item{what}{"drop" or "keep"}

  \item{pp.fold}{the minimum posterior odds for a child
  model to be returned}

  \item{quiet}{default FALSE, if TRUE, supress progress
  messages.}
}
\value{
  object of class dropModels defining models to drop (if
  what=="drop"), or an index vector of which rows of
  supplied m.child should be kept or kept (if what=="keep")
}
\description{
  mcomp.detail, internal function
}
\details{
  The prior odds are used, together with Bayes Factors, to
  determine posterior odds on which basis a subset of child
  models are returned
}
\author{
  Chris Wallace
}

