## proper data
data(for.exercise, package="snpStats")
X <- snps.10[,11:30]
n <- nrow(X)
set.seed(12346)
Y <- rnorm(n,mean=as.numeric(X[,1]))*sqrt(0.1) +
  rnorm(n,mean=as.numeric(X[,5]))*sqrt(0.1) +
  rnorm(n,mean=as.numeric(X[,11]))*sqrt(0.1) +
  rnorm(n)*sqrt(0.7)

ld <- show.ld(X=X)

(causal <- colnames(X)[c(1,5,11)])

tags <- tag(X, tag.threshold=0.8)
tags[causal]

data <- make.data(X, Y,tags=tags,family="gaussian")
bma.1 <- bma.nsnps(data, nsnps=1)
head(ss1 <- snp.summary(bma.1))
(rownames(ss1[apply(ss1,1,max) < 1, ]) -> snps.drop)

(tags <- tags[!(tags %in% snps.drop)])
data2 <- make.data(X, Y,tags=tags,family="gaussian")

bma.2 <- bma.nsnps(data, nsnps=2)


bma.2d <- bma.nsnps(data2, nsnps=2)



priors <- prior.nsnps(0:10, unique(tags), expected=3)
plot(0:10,priors, type="b")
priors <- prior.nsnps(1:10, unique(tags), expected=3)

md <- models.diff(parents=bma.1, children=bma.2d,
                  prior.parents=priors[1],
                  prior.children=priors[2])
mp <- models.prune(parents=bma.1, children=bma.2d,
                  prior.parents=priors[1],
                  prior.children=priors[2])

models.drop <- mcomp(parents=bma.1, children=bma.2d,
                     prior.parents=priors[1],
                     prior.children=priors[2])
bma.3 <- bma.nsnps(data, nsnps=3)
bma.3d <- bma.nsnps(data2, nsnps=3)
bma.3dd <- bma.nsnps(data2, nsnps=3, models.drop=models.drop)

## now check which models were dropped, plot network, was it sensible?

bma <- list(bma.1, bma.2, bma.3, bma.4)
bma <- list(bma.1, bma.2d, bma.3dd)

gplot <- function(bma) {
   priors2 <- priors[1:3]/sapply(as.list(1:3), function(i) { bma[[i]]@nmodels })
   graphs <- lapply(as.list(1:2), function(i) {
     relate(parents=bma[[i]], children=bma[[i+1]],
            name.parent=paste("M",i,sep=""), name.child=paste("M",i+1,sep=""),
            prior.parent=priors[i]/bma[[i]]@nmodels,
            prior.child=priors[i+1]/bma[[i+1]]@nmodels)
   })
  g <- graphs.merge(graphs)
   g <- add.attributes(g)
   
   df <- get.data.frame(g, what="both")
   v <- df[[1]]
   v$logpp <- v$BF + 2*log(priors2[v$x])
   v$pp <- exp(v$logpp/2)/sum(exp(v$logpp/2))
   
   from <- merge(df[[2]], v[,c("name","x","pp","logpp")], by.y="name", by.x="from")
   to <- merge(df[[2]], v[,c("name","x","pp","logpp")], by.y="name", by.x="to")
   edges <- merge(from, to, by=c("from", "to"), suffixes=c(".from",".to"))
   
   logpp <- ggplot(v, aes(x=x, y=logpp))  +
      geom_segment(data=edges, mapping=aes(x=x.from, xend=x.to, y=logpp.from, yend=logpp.to)) +
        geom_point(aes(col=logpp), size=5)
   
   pp <- ggplot(v, aes(x=x, y=pp))  +
     geom_segment(data=edges, mapping=aes(x=x.from, xend=x.to, y=pp.from, yend=pp.to)) +
       geom_point(aes(col=pp), size=5)
   return(list(logpp=logpp, pp=pp))
 }

aplots <- gplot(list(bma.1, bma.2, bma.3))
dplots <- gplot(list(bma.1, bma.2d, bma.3d))
ddplots <- gplot(list(bma.1, bma.2d, bma.3dd))

library(gridExtra)
do.call(grid.arrange,c(aplots,dplots,ddplots))
