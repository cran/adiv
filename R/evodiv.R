evodiv <-
function(phyl, comm, method = "full", tol = 1e-8){
 
  if(any(!method%in%c("richness", "GiniSimpson", "Simpson", "Shannon", "Margalef", "Menhinick", "McIntosh", "full"))) stop("Your choice for method is not available")
  if("full"%in%method) method <- c("richness", "GiniSimpson", "Simpson", "Shannon", "Margalef", "Menhinick", "McIntosh")

  tre <- .checkphyloarg(phyl)
  tre4 <- tre$phyl
  m <- comm
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
  if(any(method%in%c("GiniSimpson", "Simpson", "Shannon")))
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
      return(rowSums(m[, namestips, drop = FALSE]))
  }
  abundancesnodes <- cbind.data.frame(lapply(des, fun))
  mBabtot <- cbind(abundancesnodes, m)
  colnames(mBabtot) <- c(nodeLabels(tre4), colnames(m))
  
  branchlengths <- getEdge(tre4, colnames(mBabtot), missing = "OK")

  branchlengths <- edgeLength(tre4)[branchlengths]

  if(any(is.na(branchlengths))) stop("the lengths of some branches are missing in the phylogenetic tree; note that lengths of zero are allowed")

  tab1 <- mBabtot
  tab2 <- (t(t(mBabtot)*branchlengths))
 
  if(ncom == 1)
  composition <- as.data.frame(sweep(tab1, 1, sum(tab2), "/"))
  else
  composition <- as.data.frame(sweep(tab1, 1, rowSums(tab2), "/"))

   FUNshannon <- function(v){
    if(length(v[v>0])==1 && abs(v-1) < tol) return(0)
    else{
    w <- v[v>0]
    b <- branchlengths[v>0]
    return(as.vector(-sum(b*w*log(w))))
    }
   }
   FUNmargalef <- function(x){
    if(abs(sum(branchlengths[x>0])-1)<tol) return(0)
    else     
 return((sum(branchlengths[x>0])-1)/log(sum(branchlengths*x)))
   }
   FUNmcIntosh <- function(x){
    if(abs(sum(branchlengths*x)-sqrt(sum(branchlengths*x^2)))<tol) return(0)
    else     
 return((sum(branchlengths*x)-sqrt(sum(branchlengths*x^2)))/(sum(branchlengths*x)-sqrt(sum(branchlengths*x))))
   }

   RES <- matrix(0, nrow(composition), length(method))
   rownames(RES) <- rownames(composition)
   colnames(RES) <- method
   for(i in 1:length(method)){
      if(method[i]=="richness")
      RES[,i] <- apply(composition, 1, function(x) sum(branchlengths[x>0]))
      if(method[i]=="GiniSimpson")
      RES[,i] <- apply(composition, 1, function(x) 1-sum(branchlengths*x^2))
      if(method[i]=="Simpson")
      RES[,i] <- apply(composition, 1, function(x) 1/sum(branchlengths*x^2))
      if(method[i]=="Shannon")
      RES[,i] <- apply(composition, 1, FUNshannon)
      if(method[i]=="Margalef")
      RES[,i] <- apply(tab1, 1, FUNmargalef)
      if(method[i]=="Menhinick")
      RES[,i] <- apply(tab1, 1, function(x) sum(branchlengths[x>0])/sqrt(sum(branchlengths*x)))
      if(method[i]=="McIntosh")
      RES[,i] <- apply(tab1, 1, FUNmcIntosh)
   }
   return(RES)

}
