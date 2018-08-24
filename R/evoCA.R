evoCA <-
function(phyl, comm, scannf = TRUE, nf = 2, abundance =
TRUE){

  tre <- .checkphyloarg(phyl)
  tre4 <- tre$phyl
  m <- comm
  nsp <- ncol(m)
  if(!is.logical(abundance)) stop("Incorrect definition for parameter named abundance")
  if(is.null(colnames(m))) stop("m must have names for column")
  if(any(!colnames(m) %in%tipLabels(tre4))) stop("m contains tip names that are not available in phyl")
  if(any(m<0)) stop("m should contain nonnegative values")
  if(any(rowSums(m)==0)) stop("empty communities should be discarded")
  if(!hasEdgeLength(tre4)){
      treape <- as(tre4, "phylo")
      tre4 <- as(compute.brlen(treape, 1), "phylo4")
  }
  if(!isRooted(tre4)){
        treape <- as(tre4, "phylo")
        treape$root.edge <- 0
        tre4 <- as(treape, "phylo4")
  }

  if(!hasNodeLabels(tre4)) nodeLabels(tre4) <- names(nodeLabels(tre4))
  else{
     e <- nodeLabels(tre4)
     e[is.na(e)] <- names(e[is.na(e)])
     nodeLabels(tre4) <- e
  }

  a <- edgeLength(tre4)
  b <- a[getEdge(tre4, rootNode(tre4))]
  if(is.na(b)){
  ab <- a
  ab[getEdge(tre4, rootNode(tre4))] <- 0
  edgeLength(tre4) <- ab
  }
  tre4 <- subset(tre4, tips.exclude=tipLabels(tre4)[!tipLabels(tre4)%in%colnames(m)])

  des <- lapply(as.vector(nodeLabels(tre4)), function(x) names(descendants(tre4, x, type="tips")))
  des <- lapply(des, function(x) x[x%in%colnames(m)])
  fun <- function(namestips){
      return(rowSums(m[, namestips]))
  }
  abundancesnodes <- cbind.data.frame(lapply(des, fun))
  mBabtot <- cbind(abundancesnodes, m)
  colnames(mBabtot) <- c(nodeLabels(tre4), colnames(m))
  if(!abundance) mBabtot[mBabtot>0] <- 1

  branchlengths <- getEdge(tre4, colnames(mBabtot), missing = "OK")

  branchlengths <- edgeLength(tre4)[branchlengths]

  if(any(is.na(branchlengths))) stop("the lengths of some branches are missing in the phylogenetic tree; note that lengths of zero are allowed")

  poidsli <- rowSums(t(t(mBabtot)*branchlengths))/sum(t(t(mBabtot)*branchlengths))
  poidsco <- colSums(as.data.frame(mBabtot))/sum(t(t(mBabtot)*branchlengths))*branchlengths
  poidscoU <- colSums(as.data.frame(mBabtot))/sum(t(t(mBabtot)*branchlengths))
  Z <- diag(1/poidsli)%*%(as.matrix(mBabtot)/sum(t(t(mBabtot)*branchlengths)))%*%diag(1/poidscoU)-1
  rownames(Z) <- rownames(m)
  colnames(Z) <- colnames(mBabtot)
  X <- as.dudi(as.data.frame(Z), poidsco, poidsli, scannf = scannf, nf = nf, call = match.call(), type = "evoCA")
  attributes(X)$phy <- as(tre4, "phylo")
  return(X)

}
