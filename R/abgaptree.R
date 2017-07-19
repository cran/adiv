abgaptree <-
function(phyl, comm, exponent = 2, wcom = c("even", "speciesab"), tol = 1e-8){

    if(is.vector(comm)) stop("comm must be a numeric matrix or data frame")
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

    if(nrow(comm) < 2)
      stop("comm must have at least 2 rows")
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
  
    tabintra <- aptree(phyl.phylo, comm, exponent = exponent)
    intra <- as.matrix(tabintra)%*%wcom

    freq1 <- sweep(comm, 1, rowSums(comm), "/")
    freq1 <- apply(freq1, 2, function(x) sum(x*wcom))
    tot <- aptree(phyl.phylo, freq1, exponent = exponent, tol = tol)

    beta <- tot - intra

    restab <- data.frame(alpha = intra, beta = beta[, 1], gamma = tot[, 1])
    return(restab)

}
