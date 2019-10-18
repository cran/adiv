evouniparam <- function(phyl, comm, method = c("hill", "tsallis", "renyi"), 
    q = 2, tol = 1e-8){
 
    ow <- options("warn")
    tre <- .checkphyloarg(phyl)
    phyl1 <- tre$phyl.phylo
    A <- write.tree(phyl1)
    if(!substr(A, nchar(A)-1, nchar(A))==");"){
       phyl1 <- read.tree(text=paste0("(", substr(A,1,nchar(A)-1), ");"))
    options(warn = -1)
    }
    method <- method[1]
    reduceTree <- function(tree){
	  C <- vcv.phylo(tree, model = "Brownian")	
       if(min(diag(C))<1){
            warning("The phylogenetic tree was re-scaled so that the shortest distance from tip to root is equal to 1")
	      tree$edge.length <- tree$edge.length/min(diag(C))
       }
	  return(tree)
    }
    if(method!="hill")
      phyl1 <- reduceTree(phyl1)

    hi <- diag(vcv.phylo(phyl1, model = "Brownian"))
    phylstar <- starTree(phyl1$tip.label, hi)
    if(length(q)==1){
      vnum <- evodivparam(phyl1, comm, method, q, tol)
      vden <- evodivparam(phylstar, comm, method, q, tol)
      v <- vnum/vden
      v[vnum < tol] <- 0
      class(v) <- "evouniparam"
      options(ow)
      return(v)
    }
    if ( length(q) > 1){
      tabnum <- evodivparam(phyl1, comm, method, q, tol)$div
      tabden <- evodivparam(phylstar, comm, method, q, tol)$div
      tab1 <- as.matrix(tabnum)/as.matrix(tabden)
      tab1[as.matrix(tabnum) < tol] <- 0
      listtotale <- list()
      listtotale$q <- q
      listtotale$uni <- as.data.frame(tab1)
      class(listtotale) <- "evouniparam"
      options(ow)
      return(listtotale)
    }
}
