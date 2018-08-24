evopcachord <-
function(phyl, comm, option = c("centred", "decentred"), w = c("evoab", "even", "speciesab"), scannf = TRUE, nf = 2, abundance = TRUE){

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
      return(rowSums(m[, namestips, drop=FALSE]))
  }
  abundancesnodes <- cbind.data.frame(lapply(des, fun))
  mBabtot <- cbind(abundancesnodes, m)
  colnames(mBabtot) <- c(nodeLabels(tre4), colnames(m))

  if(!abundance) mBabtot[mBabtot>0] <- 1
  branchlengths <- getEdge(tre4, colnames(mBabtot), missing = "OK")

  branchlengths <- edgeLength(tre4)[branchlengths]

  if(any(is.na(branchlengths))) stop("the lengths of some branches are missing in the phylogenetic tree; note that lengths of zero are allowed")

  tab1 <- mBabtot
  tab2 <- (t(t(mBabtot^2)*branchlengths))

  composition <- as.data.frame(sweep(tab1, 1, sqrt(rowSums(tab2)), "/"))

  nsites <- nrow(m)
  if(is.numeric(w)){
      if(length(w)!=nsites) stop("Incorrect w")
      if(any(w<0)) stop("Positive w required")
      w <- w/sum(w)
  }
  else{
      w <- w[1]
      if(!w%in%c("evoab", "even", "speciesab")) stop("Incorrect w")
      else if(w=="evoab")   w <- rowSums(tab2)/sum(tab2)
      else if(w=="speciesab") w <- rowSums(m)/sum(m)
      else if(w=="even") w <- rep(1/nsites, nsites)
  }

  poidsco <- branchlengths
  names(poidsco) <- colnames(mBabtot)

  Z <- as.data.frame(composition)

  if(any(poidsco<1e-12)){
    Z <- Z[, poidsco>1e-12]
    tab1 <- tab1[, poidsco>1e-12]
    branchlengths <- branchlengths[poidsco>1e-12]
    poidsco <- poidsco[poidsco>1e-12]
  }

  option <- option[1]
  if(option=="centred"){
    X <- dudi.pca(Z, col.w = poidsco, row.w = w, center = TRUE, scale = FALSE, scannf = scannf, nf = nf)
  }
  else if(option=="decentred"){
    tabw1 <- t(t(tab1)*w)
    vcenter <- colSums(tabw1)/sqrt(sum(colSums(tabw1)^2*branchlengths))
    X <- dudi.pca(Z, col.w = poidsco, row.w = w, center = vcenter, scale = FALSE, scannf = scannf, nf = nf)
  }
  X$call <- match.call()
  class(X) <- c("evopca", class(X))
  attributes(X)$phy <- as(tre4, "phylo")
  return(X)

}
