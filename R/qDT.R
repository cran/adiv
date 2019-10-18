qDT <-
function(phyl, comm, q = 2, tol = 1e-8){
 
  tre <- .checkphyloarg(phyl)
  tre4 <- tre$phyl
  m <- comm
  if(!(is.numeric(q) | is.integer(q))) stop("Incorrect definition for q")
   if(any(q < -tol)) stop("q must be nonnegative")
  if(is.null(colnames(m))) stop("m must have names for column")
  if(any(colSums(m)==0)){
      nsp <- ncol(m)
      nspreal <-length((1:nsp)[colSums(m) > 0])
      if(nspreal>1)
      	m <- m[, colSums(m) > 0, drop = FALSE] 
  }
  ncom <- nrow(m)
  if(any(m<0)) stop("m should contain nonnegative values")
  if(any(rowSums(m)==0)) stop("empty communities should be discarded")
  if(any(!colnames(m)%in%tipLabels(tre4))) stop("m contains tip names that are not available in phyl")
  
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
  
  branchlengths <- getEdge(tre4, colnames(mBabtot), missing = "OK")

  branchlengths <- edgeLength(tre4)[branchlengths]

  if(any(is.na(branchlengths))) stop("the lengths of some branches are missing in the phylogenetic tree; note that lengths of zero are allowed")

  tab1 <- mBabtot
  tab2 <- (t(t(mBabtot)*branchlengths))
 
  composition <- as.data.frame(sweep(tab1, 1, rowSums(tab2), "/"))

  vresD <- rowSums(tab2)/rowSums(comm)

  hill <- function(x, q){
        funhill <- function(y, q) {
            b <- branchlengths[y>0]
            y <- y[y>0]
            if(abs(q-1) < tol){
                resi <- exp(-sum(b*y*log(y)))
            }
            else{
            	resi <- (sum(b*y^q))^(1/(1-q))
            }
            return(resi)
        }
        res <- sapply(as.data.frame(t(x)), funhill, q)
        return(res)
  }
  funq <- function(q){
        vres <- hill(composition, q)/vresD
        return(vres)
  }
  if(length(q)==1){
        v <- funq(q)
        class(v) <- "evodivparam"
        return(v)
  }
  if (length(q) > 1){
        calcul1 <- sapply(q, funq)
        if(ncom>1)
            tab1 <- cbind.data.frame(calcul1)
        else{
            tab1 <- as.data.frame(matrix(calcul1, 1, byrow=TRUE))
            rownames(tab1) <- rownames(m)
        }
        listtotale <- list()
        listtotale$q <- q
        listtotale$div <- tab1
        class(listtotale) <- "evodivparam"
        return(listtotale)
    }
}
