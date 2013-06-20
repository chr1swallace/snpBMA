graphBMA <- function(bma.list, priors) {
  n <- length(bma.list)
  priors2 <- priors[1:n]
  ## if(!is.null(prior.prob)) {
  ##   priors2 <- priors[1:n]/sapply(as.list(1:n), function(i) { bma.list[[i]]@nmodels })
  ## } else {
  ##   priors2 <- prior.odds
  ## }
  graphs <- lapply(as.list(1:(n-1)), function(i) {
    ## nparents <- if(length(bma.list[[i]]@groups)) { length(bma.list[[i]]@groups) } else { bma.list[[i]]@nmodels }
    ## nchild <- if(length(bma.list[[i+1]]@groups)) { length(bma.list[[i+1]]@groups) } else { bma.list[[i+1]]@nmodels }
    relate(parents=bma.list[[i]], children=bma.list[[i+1]],
           name.parent=paste("M",i,sep=""), name.child=paste("M",i+1,sep=""),
           pp.parent=post.snpBMA(bma.list[[i]], priors2[i]),
           pp.child=post.snpBMA(bma.list[[i+1]], priors2[i+1]),
           groups=bma.list[[i]]@groups, do.edges=1)
  })
  g <- graphs.merge(graphs)
  g <- add.attributes(g)
  return(g)
}

graphView <- function(g, what="pp") {   
   df <- get.data.frame(g, what="both")
   v <- df[[1]]
   ## v$logpp <- v$BF + 2*log(priors2[v$x])
   v$pp <- v$pp/sum(v$pp)
   v$y <- v[,what]
   
   from <- merge(df[[2]], v[,c("name","x","pp")], by.y="name", by.x="from")
   to <- merge(df[[2]], v[,c("name","x","pp")], by.y="name", by.x="to")
   edges <- merge(from, to, by=c("from", "to"), suffixes=c(".from",".to"))
   
   edges$yfrom <- edges[,paste(what, "from", sep=".")]
   edges$yto <- edges[,paste(what, "to", sep=".")]
   edges$col <- log(edges$yfrom)
   
   ggplot(v, aes(x=x, y=y))  +
##      geom_segment(data=edges, mapping=aes(x=x.from, xend=x.to, y=pp.from, yend=pp.to, linetype=type, alpha=0.05)) +
      geom_segment(data=edges, mapping=aes(x=x.from, xend=x.to, y=yfrom, yend=yto), col="grey") +
        geom_point(aes(col=log(y)), size=5) +
        geom_text(mapping=aes(label=name), hjust=0) +
          xlim(c(min(v$x), max(v$x) + 0.3))
}

within.graph <- function(models, groups, name, target) {
  mg.models <- models.group.collapse(models, groups=groups)
  models.linked <- crossprod(t(mg.models))
  tmp <- Matrix(0,nrow(models.linked),ncol(models.linked),
                dimnames=list(paste(name,model.names(mg.models),sep="-"),
                  paste(name,model.names(mg.models),sep="-")))
  tmp[models.linked==target] <- 1
  diag(tmp) <- 0
  g.child <- graph.adjacency(tmp)
  g.child <- set.edge.attribute(g.child, "LD", value=TRUE)
  return(g.child)
}

model.names <- function(models) {
  m <- as(models, "dgTMatrix")
  sapply(split(m@j,m@i), paste, collapse="-")
}

relate <- function(parents, children, name.parent="parent", name.child="child",
                   pp.parent=NULL, pp.child=NULL, groups=list(), do.edges=0) {
  m.parent=parents@models
  m.child=children@models
  bf.parent=parents@bf
  bf.child=children@bf
  ntarget=parents@nsnps
  cols.use <- intersect(colnames(m.parent),colnames(m.child))
  relate <- m.parent[,cols.use] %*% Matrix::t(m.child[,cols.use])
  index <- Matrix::which(relate==ntarget, arr.ind=TRUE)
#  rel.prior <- 2*(log(prior.children) - log(prior.parents))
  relate2 <- Matrix(0,nrow(relate),ncol(relate),sparse=TRUE)
  relate2[index] <- bf.parent[ index[,1], 2 ] - bf.child[ index[,2], 2]
  
  dimnames(relate2) <- list(paste(name.parent,model.names(m.parent),sep="-"),
                            paste(name.child,model.names(m.child),sep="-"))
  g <- graph.incidence(relate2)

  ## group edges
  ## if(length(groups)) {
  ##    gList <- list(set.edge.attribute(g, "LD", value=FALSE))
  ##   if(do.edges==2) ## children
  ##     gList <- c(gList, list(within.graph(m.child, groups, name.child, children@nsnps)))
  ##   if(do.edges>0) ## parents
  ##     gList <- c(gList, list(within.graph(m.parent, groups, name.parent, parents@nsnps)))
  ##   g <- graphs.merge(gList)
  ##    g <- edge.colour(g, "LD")
  ## }
                             
  g <- set.vertex.attribute(g, "generation", rownames(relate2), name.parent)
  g <- set.vertex.attribute(g, "generation", colnames(relate2), name.child)
  g <- set.vertex.attribute(g, "BF", rownames(relate2), bf.parent[,2])
  g <- set.vertex.attribute(g, "BF", colnames(relate2), bf.child[,2])
  if(!is.null(pp.parent))
    g <- set.vertex.attribute(g, "pp", rownames(relate2), pp.parent[,2])
  if(!is.null(pp.child))
    g <- set.vertex.attribute(g, "pp", colnames(relate2), pp.child[,2])
  return(g)
}

vertex.colour <- function(g, name, n=20) {
  v <- get.vertex.attribute(g, name)
  cols <- heat.colors(n)
  vc <- cut(v, n)
  set.vertex.attribute(g, "color", index=V(g), value=cols[ as.numeric(vc) ])
}

edge.colour <- function(g, name) {
  v <- get.edge.attribute(g, name)
  cols <- ifelse(v, "red", "black")
  set.edge.attribute(g, "color", index=E(g), value=cols)
}

vertex.shape <- function(g, name) {
  shape.names <- c("circle", "square", "csquare", "rectangle", "crectangle",
                   "vrectangle", "pie")  
  v <- factor(get.vertex.attribute(g, name))
  set.vertex.attribute(g, "shape", index=V(g), value=shape.names[ as.numeric(v) ])
}

vertex.size <- function(g, name, offset=5) {
  v <- factor(get.vertex.attribute(g, name))
  sizes <- seq_along( levels(v) ) + offset
  set.vertex.attribute(g, "size", index=V(g), value=sizes[ as.numeric(v) ])
}

vertex.xy <- function(g, xname, yname) {
  x <- as.numeric(factor(get.vertex.attribute(g, xname)))
  y <- get.vertex.attribute(g, yname)
  g <- set.vertex.attribute(g, "x", index=V(g), value=x)
  g <- set.vertex.attribute(g, "y", index=V(g), value=y)
}

add.attributes <- function(g) {
  g <- vertex.colour(g, "pp")
#  g <- edge.colour(g, "weight")
#  g <- vertex.size(g, "generation")
  g <- vertex.xy(g, "generation", "pp")
}

graphs.merge <- function(gList) {
  dfList <- lapply(gList, get.data.frame, what="both")
  eList <- lapply(dfList, "[[", 2)
  vList <- lapply(dfList, "[[", 1)
  vList <- lapply(vList, function(v) v[, setdiff(colnames(vList[[1]]), "type"), drop=FALSE])
  e <- do.call("rbind",eList)
  v <- unique(do.call("rbind", vList))
  g <- graph.data.frame(d=e, directed=FALSE, vertices=v)
  return(g)
}
