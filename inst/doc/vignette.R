### R code from vignette source 'vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: vignette.Rnw:47-57
###################################################
library(snpStats)
data(for.exercise, package="snpStats")
X <- snps.10[,11:30]
n <- nrow(X)
set.seed(12346)
Y <- rnorm(n,mean=as.numeric(X[,1]))*sqrt(0.1) +
  rnorm(n,mean=as.numeric(X[,5]))*sqrt(0.1) +
  rnorm(n,mean=as.numeric(X[,11]))*sqrt(0.1) +
  rnorm(n)*sqrt(0.7)
(causal <- colnames(X)[c(1,5,11)])


###################################################
### code chunk number 2: vignette.Rnw:65-67
###################################################
summary(X)
X <- impute.missing(X)


###################################################
### code chunk number 3: vignette.Rnw:74-75
###################################################
ld <- show.ld(X=X)


###################################################
### code chunk number 4: vignette.Rnw:95-97
###################################################
tags <- tag(X, tag.threshold=0.8)
tags[causal]


###################################################
### code chunk number 5: vignette.Rnw:103-111
###################################################
## make a snpBMAdata set
data <- make.data(X, Y,tags=tags,family="gaussian")

## Calculate Bayes Factors for all one SNP models
bma.1 <- bma.nsnps(data, nsnps=1)

## Summarise the SNPs with greatest support
head(ss1 <- snp.summary(bma.1))


###################################################
### code chunk number 6: vignette.Rnw:119-122
###################################################
bma.2 <- bma.grow(data=data, bma=bma.1)
bma.3 <- bma.grow(data=data, bma=bma.2)
bma.4 <- bma.grow(data=data, bma=bma.3)


###################################################
### code chunk number 7: vignette.Rnw:138-141
###################################################
## assume a binomial prior for the number of SNPs with expectation of 3 causal SNPs
## ie exactly the scenario simulated!
priors <- prior.binomial(1:10, n=ncol(X), expected=3)


###################################################
### code chunk number 8: vignette.Rnw:147-156
###################################################
## create a graph of BMA results so far
results <- stack(bma.1,bma.2,bma.3,bma.4)
 g<-graphBMA(bma.list=results, priors)

## g is an igraph, so you can do all the usual stuff with it:
g

## visualize
graphView(g)


###################################################
### code chunk number 9: vignette.Rnw:164-165
###################################################
snps0(bma.3)[ as.character(c(0,1,5)) ]


###################################################
### code chunk number 10: vignette.Rnw:169-170
###################################################
top.models(results, priors)


###################################################
### code chunk number 11: vignette.Rnw:191-198
###################################################
## First, tag at r2=0.99
tags.99 <- tag(X, 0.99)

## group remaining snps at r2=0.8, using the first set of tags above as indices
groups <- group.tags(tags, keep=tags.99)
length(groups)
data.99 <- make.data(X, Y, tags=tags.99, family="guassian")


###################################################
### code chunk number 12: vignette.Rnw:205-207
###################################################
top.models(results, priors)
expand.snps <- top.snps(results, priors, nmodels=3)


###################################################
### code chunk number 13: vignette.Rnw:212-216
###################################################
bma.e1 <- bma.expand(data.99, bma.1, groups=groups[expand.snps])
bma.e2 <- bma.expand(data.99, bma.2, groups=groups[expand.snps])
bma.e3 <- bma.expand(data.99, bma.3, groups=groups[expand.snps])
bma.e4 <- bma.expand(data.99, bma.4, groups=groups[expand.snps])


###################################################
### code chunk number 14: vignette.Rnw:222-230
###################################################
## create a graph of BMA results so far
expand.results <- stack(bma.e1,bma.e2,bma.e3,bma.e4)
 g.expand<-graphBMA(expand.results, priors)

## visualize
graphView(g.expand)

top.models(expand.results, priors)


