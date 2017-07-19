aptree <-
function(phyl, comm, exponent = 2, tol = 1e-8){
 
    if(is.vector(comm)) comm <- t(comm)
    if(any(rowSums(comm)<tol)) stop("empty communities must be removed")
    tre <- .checkphyloarg(phyl)
    phyl.phylo <- tre$phyl.phylo

    if(length(phyl.phylo$tip.label)==ncol(comm)){
        if(!is.null(colnames(comm)) & any(!phyl.phylo$tip.label%in%colnames(comm))) stop("names of species in comm are not equal to tip names in phyl")
    }
    else if(length(phyl.phylo$tip.label)<ncol(comm)) stop("phyl must contain all species in comm")
    else{
        if(any(!colnames(comm)%in%phyl.phylo$tip.label)) stop("some names of species in comm are not in tip names of phyl")
        else
            phyl.phylo <- drop.tip(phyl.phylo, phyl.phylo$tip.label[!phyl.phylo$tip.label%in%colnames(comm)])
    }    
    if(!is.null(colnames(comm)))
        comm <- comm[, phyl.phylo$tip.label, drop=FALSE]
        
    if(!all(apply(comm, 1, is.numeric))) stop("comm must be a numeric matrix or data frame")
    if(any(comm < (-tol))) stop("comm must have nonnegative values")
  
    shannon <- function(x) ifelse(sum(x)>1e-8, -sum(x[x>0]/sum(x)*log(x[x>0]/sum(x),2), 0))
    charvat <- function(x, i) (1-sum(x[x>0]^i))/(i-1)
    richnessminusone <- function(x) length(x[x>0])-1

    phylD <- cophenetic.phylo(phyl.phylo) / 2
    treeh <- hclust(as.dist(phylD), method = "average")
  
    fun <- function(freq1){
    
        freq1 <- freq1/sum(freq1)
        snods <- sort(as.vector(phylD))
        dnods <- diff(snods)
        if(tol > min(dnods)) warnings("the tolerance is too high")
        vnods <- c(snods[1], snods[-1][dnods > tol])
  
        ncut <- length(vnods)

        vnodsbis <- c(0, vnods)
        dvnodsbis <- diff(vnodsbis)
        listcut <- lapply(1:ncut, function(i) cutree(treeh, h = vnods[i] - tol))
        abcut <- lapply(listcut, function(x) tapply(freq1, as.factor(x), sum))
        if(abs(exponent - 1) < tol)
            simcut <- unlist(lapply(abcut, shannon))
        else if(exponent < tol)
            simcut <- unlist(lapply(abcut, richnessminusone))
        else
            simcut <- unlist(lapply(abcut, charvat, exponent))
    
        vnodsbis <- c(0, vnods)
        dvnodsbis <- diff(vnodsbis)
        valcut <- dvnodsbis*simcut
        return(valcut)
    }
    restab <- as.data.frame(apply(comm, 1, fun))
    class(restab) <- c("aptree", "data.frame")
    attributes(restab)$call <- match.call()
    return(restab)
  
}
