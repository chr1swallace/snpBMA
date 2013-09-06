\name{r2.groups}
\alias{r2.groups}
\title{create groups of SNPs with LD > r2.threshold with index.snps}
\usage{
  r2.groups(X, r2.threshold = 0.9, snps = NULL,
    samples = NULL, do.plot = FALSE)
}
\arguments{
  \item{X}{SnpMatrix}

  \item{r2.threshold}{r2 threshold above which SNPs should
  be tagged, default 0.9}

  \item{snps}{optional, if present, limit to these SNPs}

  \item{samples}{optional, if present, limit to these
  samples}

  \item{do.plot}{default FALSE, if TRUE plot the hclust
  object used for the clustering}
}
\description{
  create groups of SNPs with LD > r2.threshold with
  index.snps
}

