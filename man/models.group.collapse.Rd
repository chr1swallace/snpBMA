\name{models.group.collapse}
\alias{models.group.collapse}
\title{models.group.collapse}
\usage{
  models.group.collapse(models, groups)
}
\arguments{
  \item{models}{}

  \item{groups}{}
}
\value{
  A dgCMatrix, with ncol==length(groups), and entry 1 if a
  SNP from that group is in the model.  Returns the input
  matrix in the case groups has zero length.
}
\description{
  Collapse a model matrix by SNP groups
}
\details{
  .. content for \details{} ..
}
\author{
  Chris Wallace
}

