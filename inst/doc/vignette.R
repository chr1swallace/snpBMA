### R code from vignette source 'vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: prior1
###################################################
library(snpBMA)
## binomial prior assuming 100 SNPs in a region, 3 expected in preferred model
## low-ish prior (0.1) of no association in region
p.binom <- prior.binomial(x=0:10, n=100, expected=3, value="prob", pi0=0.1)
## beta binomial, as above, allowing some overdispersion
p.bb.1.1 <- prior.betabinomial(x=0:10, n=100, expected=3, 
                               overdispersion=1.1, value="prob", pi0=0.1)
p.bb.1.5 <- prior.betabinomial(x=0:10, n=100, expected=3, 
                               overdispersion=1.5, value="prob", pi0=0.1)

## plot
plot(0:10, p.binom, type="b", main="Prior for the number of SNPs in chosen model",
     xlab="Number of SNPs", ylab="Prior probability")
lines(0:10, p.bb.1.1, col="green", type="b")
lines(0:10, p.bb.1.5, col="blue", type="b")
legend("topright",lty=rep(1,3), 
       col=c("black","green","blue"), 
       legend=c("Binomial", "BBinom, 1.1", "BBinom, 1.5"))


###################################################
### code chunk number 2: vignette.Rnw:91-94
###################################################
mx <- sapply(0:10, function(i) max.models(n.snps=100, n.use=i))
plot(0:10, log(mx), type="b", ylab="Log number of models",
     main="Number of models for a fixed number of SNPs (log scale)",xlab="Number of SNPs")


###################################################
### code chunk number 3: fig
###################################################
plot(0:10, p.bb.1.5/mx, type="b", main="Prior for each model of a given size",col="blue",
     xlab="Number of SNPs", ylab="Prior probability")
lines(0:10, p.bb.1.1/mx, col="green", type="b")
lines(0:10, p.binom/mx, col="black", type="b")
legend("topright",lty=rep(1,3), 
       col=c("black","green","blue"), 
       legend=c("Binomial", "BBinom, 1.1", "BBinom, 1.5"))


###################################################
### code chunk number 4: fig
###################################################
## binomial prior assuming 100 SNPs in a region, 3 expected in preferred model
odds.binom <- prior.binomial(x=0:10, n=100, expected=3, value="odds", pi0=0.1)
## beta binomial, as above, allowing some overdispersion
odds.bb.1.1 <- prior.betabinomial(x=0:10, n=100, expected=3, 
                               overdispersion=1.1, value="odds", pi0=0.1)
odds.bb.1.5 <- prior.betabinomial(x=0:10, n=100, expected=3, 
                               overdispersion=1.5, value="odds", pi0=0.1)


plot(0:10, odds.binom, type="b", main="Prior odds for each model of a given size")
lines(0:10, odds.bb.1.1, col="green", type="b")
lines(0:10, odds.bb.1.5, col="blue", type="b")
legend("topright",lty=rep(1,3), 
       col=c("black","green","blue"), 
       legend=c("Binomial", "BBinom, 1.1", "BBinom, 1.5"))


###################################################
### code chunk number 5: sim
###################################################
library(snpStats)
data(for.exercise, package="snpStats")
set.seed(12346)
X <- snps.10[,101:120]
n <- nrow(X)
causal <- c("rs1555897","rs7069505")
Y <- rnorm(n,mean=as.numeric(X[,causal[1]]))*sqrt(0.2) +
  rnorm(n,mean=as.numeric(X[,causal[2]]))*sqrt(0.2) +
  rnorm(n)*sqrt(0.6)


###################################################
### code chunk number 6: vignette.Rnw:155-158
###################################################
summary(col.summary(X))
X <- impute.missing(X)
summary(col.summary(X))


###################################################
### code chunk number 7: vignette.Rnw:165-166
###################################################
ld <- show.ld(X=X)


###################################################
### code chunk number 8: tag
###################################################
tags <- tag(X, tag.threshold=0.8)
tags[causal]


###################################################
### code chunk number 9: bma1
###################################################
## make a snpBMAdata set
data <- make.data(X, Y,tags=tags,family="gaussian")

## Calculate Bayes Factors for all one SNP models
bma.1 <- bma.nsnps(data, nsnps=1)

## Summarise the SNPs with greatest support
head(ss1 <- snp.summary(bma.1))


###################################################
### code chunk number 10: bma2
###################################################
bma.2 <- bma.grow(data=data, bma=bma.1)
bma.3 <- bma.grow(data=data, bma=bma.2)


###################################################
### code chunk number 11: vignette.Rnw:230-236
###################################################
  ## assume a binomial prior for the number of SNPs with expectation of
  ## 2 causal SNPs ie exactly the scenario simulated!  
prior.odds <- prior.binomial(1:10, n=ncol(X), expected=2, value="odds", pi0=0.1)
  ## also consider a prior with expectation of 4 causal
  ## SNPs to examine the effect of varying the prior
prior.alt <- prior.binomial(1:10, n=ncol(X), expected=4, value="odds", pi0=0.1)


###################################################
### code chunk number 12: vignette.Rnw:242-251
###################################################
## create a graph of BMA results so far
results <- stack(bma.1,bma.2,bma.3)
 g<-graphBMA(bma.list=results, priors=prior.odds)

## g is an igraph, so you can do all the usual stuff with it:
g

## visualize
graphView(g)


###################################################
### code chunk number 13: top.models
###################################################
top.models(results, priors=prior.odds, n=20)


###################################################
### code chunk number 14: vignette.Rnw:284-291
###################################################
## First, tag at r2=0.99
tags.99 <- tag(X, 0.99)

## group remaining snps at r2=0.8, using the first set of tags above as indices
groups <- group.tags(tags, keep=tags.99)
length(groups)
data.99 <- make.data(X, Y, tags=tags.99, family="guassian")


###################################################
### code chunk number 15: vignette.Rnw:300-304
###################################################
expand.snps <- top.snps(results, prior.odds, nmodels=17)
bma.e1 <- bma.expand(data.99, bma.1, groups=groups[expand.snps])
bma.e2 <- bma.expand(data.99, bma.2, groups=groups[expand.snps])
bma.e3 <- bma.expand(data.99, bma.3, groups=groups[expand.snps])


###################################################
### code chunk number 16: vignette.Rnw:312-320
###################################################
## create a graph of BMA results so far
expand.results <- stack(bma.e1,bma.e2,bma.e3)
 g.expand<-graphBMA(expand.results, prior.odds)

## visualize
graphView(g.expand)

top.models(expand.results, prior.odds)


###################################################
### code chunk number 17: vignette.Rnw:331-335
###################################################
## define the list of SNPs to drop
max.bf <- apply(ss1,1,max)
snps.drop <- rownames(ss1)[ max.bf < 0 ]
snps.drop


###################################################
### code chunk number 18: vignette.Rnw:343-352
###################################################
## generate a new set of tags and snpBMAdata object
data2 <- snps.prune.data(data, snps.drop)

bma.2 <- bma.nsnps(data, nsnps=2)
bma.2d <- bma.nsnps(data2, nsnps=2)

## compare top models
top.models(bma.2)
top.models(bma.2d)


###################################################
### code chunk number 19: BMA3
###################################################
## prune the bma.2d object
bma.2dd <- models.prune(parents=bma.1, children=bma.2d, 
                        prior.parents=prior.odds[1],
                        prior.children=prior.odds[2])

## grow the BMA to a third generation
bma.3dd <- bma.grow(data2, bma.2dd)

## for comparison, without pruning, we could use tagging only...
bma.3 <- bma.nsnps(data, nsnps=3)

## ... or tagging + excluding poorly supported single SNPs
bma.3d <- bma.nsnps(data2, nsnps=3)

## this should be the same as growing from the bma.2d object
bma.3d2 <- bma.grow(data2, bma.2d)

top.models(bma.3d)
top.models(bma.3d2)


###################################################
### code chunk number 20: vignette.Rnw:433-435
###################################################
sdata <- make.data(X,Y,tags=tags,family="gaussian",
                  strata=rep(c(1,2),length=nrow(X)))


###################################################
### code chunk number 21: vignette.Rnw:442-448
###################################################
data <- make.data(X,Y,tags=tags,family="gaussian")
bma.1 <- bma.nsnps(data, nsnps=1)
sbma.1 <- bma.nsnps(sdata, nsnps=1)

top.models(bma.1)
top.models(sbma.1)


