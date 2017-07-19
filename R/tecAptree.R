tecAptree <-
function(phyl, v = NULL, tol = 1e-8){
  
    tre <- .checkphyloarg(phyl)
    phyl.phylo <- tre$phyl.phylo
    
    if(!is.null(v)){
        if(!is.vector(v)) stop("v must be a vector")
        if(!is.numeric(v)) stop("v must be numeric")
        if(length(phyl.phylo$tip.label) != length(v))
            stop("the number of values in v must be equal to the number of tips in phyl")
    }

    phyl.D <- cophenetic.phylo(phyl.phylo) / 2
    snods <- sort(as.vector(phyl.D))
    dnods <- diff(snods)
    if(tol > min(dnods)) warnings("the tolerance, tol, is too high")
    vnods <- c(snods[1], snods[-1][dnods > tol])
  
    ncut <- length(vnods)
    treeh <- hclust(as.dist(phyl.D), method = "average")
    vnodsbis <- c(0, vnods)
    dvnodsbis <- diff(vnodsbis)
    listcut <- lapply(1:ncut, function(i) cutree(treeh, h = vnods[i] - tol))

    names(listcut) <- paste("p", 1:length(listcut), sep = "")
    riccut <- unlist(lapply(listcut, function(x) length(unique(x))))
    if(!is.null(v)){
        abcut <- lapply(listcut, function(x) tapply(v, as.factor(x), sum)) 
        names(abcut) <- paste("p", 1:length(abcut), sep = "")
    }
    vnodsbis <- c(0, vnods)
    dvnodsbis <- diff(vnodsbis)
    if(!is.null(v))
        res <- list(h = vnods, plength = dvnodsbis, nbgroups = riccut, list = listcut, relab = abcut, call = match.call())
    else
        res <- list(h = vnods, plength = dvnodsbis, nbgroups = riccut, list = listcut, call = match.call())
    return(res)

}
