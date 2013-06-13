\name{index.groups}
\alias{index.groups}
\title{index.groups}
\usage{
  index.groups(X, index.snps, r2.threshold = 0.9,
    snps = NULL, samples = NULL)
}
\arguments{
  \item{X}{SnpMatrix object}

  \item{index.snps}{snps to index groups}

  \item{r2.threshold}{threshold above which snps should be
  grouped with an index snp}

  \item{snps}{SNP subset to be used (columns of X)}

  \item{samples}{sample subset to be used (rows of X)}
}
\value{
  list of character vectors, each vector being a single
  group and the names of the list being the index snps
}
\description{
  create groups of SNPs with LD > r2.threshold with
  index.snps
}
\details{
  .. content for \details{} ..
}
\author{
  Chris Wallace
}

