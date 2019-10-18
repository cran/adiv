evoeveparam <-
function(phyl, comm, method = c("hill", "tsallis", "renyi"), 
    q = 2, option = 1:3, H = NULL, tol = 1e-8){
 
  ow <- options("warn")
  tre <- .checkphyloarg(phyl)
  PHYTEST <- tre$phyl.phylo
  A <- write.tree(PHYTEST)
  if(!substr(A, nchar(A)-1, nchar(A))==");"){
     PHYTEST <- read.tree(text=paste0("(", substr(A,1,nchar(A)-1), ");"))
  options(warn = -1)
  }
  tre4 <- as(PHYTEST, "phylo4")
  m <- comm
  method <- method[1]
  option <- option[1]
  if(!option%in%(1:3)) stop("Unavailable option")
  if(!method%in%c("tsallis", "hill", "renyi")) stop("Unavailable method")
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
  
  reduceTree <- function(tree){
	  C <- vcv.phylo(tree, model = "Brownian")	
       if(min(diag(C))<1){
            warning("The phylogenetic tree was re-scaled so that the shortest distance from tip to root is equal to 1")
	      tree$edge.length <- tree$edge.length/min(diag(C))
       }
	  return(tree)
  }

  if(method!="hill")
      tre4 <- as(reduceTree(as(tre4, "phylo")), "phylo4")
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
  tre.phylo <- as(tre4, "phylo")
  options(ow)

  V <- diag(vcv.phylo(tre.phylo, model = "Brownian"))
  names(V) <- tre.phylo$tip.label
  mV <- m[, names(V)]
  if(option==1)
      CC <- apply(mV, 1, function(x) sum(V[x>0]))
  if(option==2)
      CC <- apply(mV, 1, function(x) length(x[x>0])*max(V[x>0]))
  if(option==3){
      if(is.null(H)){
         warning("Parameter H, a numeric, must be specified if option = 3; option 2 was used here")
         CC <- apply(mV, 1, function(x) length(x[x>0])*max(V[x>0]))
      }
      else
      CC <- apply(mV, 1, function(x) length(x[x>0])*H)
  }

  tsallis <- function(x, q){
        funtsallis <- function(y, q) {
            b <- branchlengths[y>0]
            y <- y[y>0]
            if(abs(q-1) < tol){
                resi <- -sum(b*y*log(y))
            }
            else{
            	resi <- (1-sum(b*y^q))/(q-1)
            }
            return(resi)
        }
        res <- sapply(as.data.frame(t(x)), funtsallis, q)
        return(res)
   }
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
    renyi <- function(x, q){
        funrenyi <- function(y, q) {
            b <- branchlengths[y>0]
            y <- y[y>0]
            if(abs(q-1) < tol){
                resi <- -sum(b*y*log(y))
            }
            else{
            	resi <- log((sum(b*y^q))^(1/(1-q)))
            }
            return(resi)
        }
        res <- sapply(as.data.frame(t(x)), funrenyi, q)
        return(res)
    }
    funq <- function(q){
        if(method == "tsallis"){
            T <- tsallis(composition, q)
            if(abs(q-1) < tol){
                vres <- T/log(CC)
            }
            else{
                vres <- T/(1-CC^(1-q))*(q-1)
            }
            vres[T < tol] <- 0
            return(vres)
        }
        if(method == "hill"){
            vres <- hill(composition, q)/CC
            return(vres)
        }
        if(method == "renyi"){
            R <- renyi(composition, q)
            vres <- R/log(CC)
            vres[R < tol] <- 0
            return(vres)
        }
    }
    if(length(q)==1){
        v <- funq(q)
        class(v) <- "evoeveparam"
        return(v)
    }
    if ( length(q) > 1){
           calcul1 <- sapply(q, funq)
           if(ncom>1)
               tab1 <- cbind.data.frame(calcul1)
           else{
               tab1 <- as.data.frame(matrix(calcul1, 1, byrow=TRUE))
               rownames(tab1) <- rownames(m)
           }
           listtotale <- list()
           listtotale$q <- q
           listtotale$eve <- tab1
           class(listtotale) <- "evoeveparam"
           return(listtotale)
    }
}
