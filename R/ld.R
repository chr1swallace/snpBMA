##' Calculate D' and R^2 for a SnpMatrix object and display as a heatmap
##'
##' Side effect: displays plot on current graphics device
##' @title show.ld
##' @return invisibly returns matrix with D' in upper.tri() entries and R^2 in lower.tri() entries
##' @author Chris Wallace
##' @export
##' @param X a SnpMatrix object
##' @param snps optional character vector of column names of X (SNPs)
##' for which LD should be calculated
##' @param samples optional character vector of row names of X
##' (Samples) for which LD should be calculated
##' @param lines.limit draw faint lines to make separation between SNPs clearer if the number of SNPs is < lines.limit
show.ld <- function(X, snps=colnames(X), samples=rownames(X), 
                    lines.limit=20) {
 my.colours <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
                 "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  if(!is(X,"SnpMatrix"))
    stop("X must be a SnpMatrix")
 breaks <- NULL
 groups <- NULL
  if(is.list(snps)) {
    breaks <- sapply(snps,length)
    groups <- rep(1:length(snps), times=breaks)
    snps <- unlist(snps)
  }
  LD <- ld(X[samples,snps],
           depth=length(snps)-1,
           stat=c("D.prime","R.squared"),
           symmetric=TRUE)
  ld <- as.matrix(LD$R.squared)
 wh <- which(upper.tri(ld))
  ld[wh] <- as.matrix(LD$D.prime)[ wh ]
 
 diag(ld) <- 1
 if(!is.null(groups))
   colnames(ld) <- rownames(ld) <- paste(groups,colnames(ld),sep="/")
 
  df <- melt(as.matrix(ld))
  df$X1 <- factor(df$X1, levels=colnames(ld))
  df$X2 <- factor(df$X2, levels=colnames(ld))
 n <- max(as.numeric(df$X1))
 p <- ggplot(df) + geom_tile(mapping=aes(x=X1,y=X2,fill=value), linetype=1) +
    scale_fill_gradientn(colours = my.colours) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text = element_text(colour="black"),
            axis.title = element_blank())
 if(n<lines.limit) {
   p <- p + geom_hline(yintercept=(0:n)+0.5, col="grey60") +
     geom_vline(xintercept=(0:n)+0.5, col="grey60")
 }
 
if(!is.null(breaks)) {
  x <- c(0,cumsum(breaks),max(as.numeric(df$X1))) + 0.5
  n <- length(x)
  xf <- rbind(data.frame(x=x[-n], y=x[-n], xend=x[-1], yend=x[-n]),
              data.frame(x=x[-n], y=x[-n], xend=x[-n], yend=x[-1]),
              data.frame(x=x[-1], y=x[-1], xend=x[-1], yend=x[-n]),
              data.frame(x=x[-1], y=x[-1], xend=x[-n], yend=x[-1]))
              
  show(p + geom_segment(aes(x=x,xend=xend,y=y,yend=yend), data=xf, size=1))
 } else {
   show(p)
 }
  ## pheatmap(ld,
  ##          cluster_rows=FALSE, cluster_cols=FALSE)
  invisible(ld)
}
myr2 <- function(X) {
  r2 <- ld(X,
           depth=ncol(X)-1,
           symmetric=TRUE,
           stat="R.squared")
  if(any(is.na(r2))) {
    r2.na <- as(is.na(r2),"matrix")
    use <- rowSums(r2.na)>0
    ## work around for r2=NA bug.  
    r2.cor <- as(cor(as(X[,use,drop=FALSE],"numeric"), use="pairwise.complete.obs")^2,"Matrix")
    r2[ which(r2.na) ] <- r2.cor[ which(r2.na[use,use]) ]
  }
  diag(r2) <- 1
  return(r2)
}
##' Derive tag SNPs for a SnpMatrix object using heirarchical clustering
##'
##' Uses complete linkage and the \code{\link{hclust}} function to define clusters, then cuts the tree at 1-tag.threshold
##' @title tag
##' @param snps colnames of the SnpMatrix object to be used
##' @param tag.threshold threshold to cut tree, default=0.99
##' @param samples optional, subset of samples to use
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @author Chris Wallace
##' @export
tag <- function(X,tag.threshold=0.99, snps=NULL, samples=NULL) {
  if(!is.null(snps) || !is.null(samples))
    X <- X[samples,snps]
  r2 <- myr2(X)
   D <- as.dist(1-r2)
   hc <- hclust(D, method="complete")
   clusters <- cutree(hc, h=1-tag.threshold)
   snps.use <- names(clusters)[!duplicated(clusters)]
   r2.use <- r2[snps.use, colnames(X), drop=FALSE]
   tags <- rownames(r2.use)[apply(r2.use,2,which.max)]
   names(tags) <- colnames(r2.use)
   return(tags)
}

group.tags <- function(tags, keep) {
  groups <- tags[ names(tags) %in% keep ]
  groups <- split(names(groups), groups)
}
