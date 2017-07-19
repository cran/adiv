apd <-
function(phyl, comm, wcom = c("even", "speciesab"), nrep = 99, alter = "two-sided", tol = 1e-8){
  
    if(is.vector(comm)) comm <- t(comm)

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
  
    D <- cophenetic.phylo(phyl.phylo)/2

    if(nrow(comm)>1){    
    if(wcom[1] == "even")
        comm <- sweep(comm, 1, rowSums(comm), "/")
    else if(wcom[1] == "speciesab")
        comm <- comm
    else if(is.numeric(wcom) & length(wcom)==nrow(comm)){
        if(any(wcom < -tol))
            stop("negative values in wcom")
        wcom <- wcom/sum(wcom)
        comm <- sweep(comm, 1, rowSums(comm), "/")
        comm <- sweep(comm, 1, wcom, "*")
    }
    else stop("incorrect definition of wcom")
    }
    freq <- colSums(comm)/sum(comm)
    
    ntips <- nTips(phyl.phylo)

    stat <- function(x){
        tot <- sum(D)/ntips/(ntips-1)
        wt <- (t(x)%*%D%*%x)/(1-sum(x^2))
        return((tot-wt)/tot)
    }
    sim <- function(i, u){
        v <- u[sample(ntips)]
        return(stat(v))
    }
    theo <- sapply(1:nrep, sim, freq)
    obs <- as.vector(stat(freq))
    res <- as.randtest(sim = theo, obs = obs, alter = alter)
    res$call <- match.call()
    return(res)
  
}
