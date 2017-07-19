rtestaptree <-
function(phyl, comm, nrep = 99, alter = "two-sided", exponent = 2, 
                         wcom = c("even", "speciesab"), tol = 1e-8){
                         
    if(!inherits(comm, "data.frame") & !inherits(comm, "matrix")) stop("comm must be a numeric matrix or data frame")
    if(nrow(comm) < 2)
      stop("comm must have at least 2 rows")
      
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
        comm <- comm[, phyl.phylo$tip.label]
    if(!all(apply(comm, 1, is.numeric))) stop("comm must be a numeric matrix or data frame")
    if(any(comm < (-tol))) stop("comm must have nonnegative values")

    if(wcom[1] == "even")
        wcom <- rep(1/nrow(comm), nrow(comm))
    else if(wcom[1] == "speciesab")
        wcom <- rowSums(comm)/sum(comm)
    else if(is.numeric(wcom) & length(wcom)==nrow(comm)){
        if(any(wcom < -tol))
            stop("negative values in wcom")
        wcom <- wcom/sum(wcom)
    }
    else stop("incorrect definition of wcom")
    
    nsp <- ncol(comm)
  
    obs <- abgaptree(phyl.phylo, comm, exponent = exponent, wcom = wcom, tol = tol)[-(1:2), ]
    obs <- obs[, 2] / obs[, 3]
    #orderobs <- rev(order(obs))
  
    funsim <- function(i){
        e <- sample(1:nsp)
        comsim <- comm[, e]
        colnames(comsim) <- colnames(comm)
        theo <- abgaptree(phyl.phylo, comsim, exponent = exponent, wcom = wcom, tol = tol)[-(1:2), ]
        theo <- theo[, 2] / theo[, 3]
        return(theo)
    }
    theotab <- t(cbind.data.frame(sapply(1:nrep, funsim)))
    #nam <- paste("p", 3:(length(obs)+2), sep = "")[orderobs]
    nam <- paste("p", 3:(length(obs)+2), sep = "")
    #res <- as.krandtest(theotab[, orderobs], obs[orderobs], alter = alter, names = nam, call = match.call())
    res <- as.krandtest(theotab, obs, alter = alter, names = nam, call = match.call())
 
    class(res) <- c("rtestaptree", class(res))
    return(res)

}
