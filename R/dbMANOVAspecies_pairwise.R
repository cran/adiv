dbMANOVAspecies_pairwise <- function(dbobj, signif = TRUE, salpha = 0.05, nrep = NULL){

    tol <- dbobj$tol
    salpha <- salpha[1]
    signif <- signif[1]
    if(salpha > 1 | salpha < 0) stop("salpha must be a numeric value between 0 and 1")
    if(!inherits(signif, "logical")) stop("Incorrect definition of argument signif")
    if(!inherits(dbobj, "dbMANOVAspecies")) stop("Argument dbobj must be of class dbMANOVAspecies")
    padjust <- dbobj$test$adj.method
    if(is.null(nrep)) nrep <- max(dbobj$test$rep)
    LAAA <- as.list(dbobj$call)
    comm <- eval.parent(LAAA$comm)
    sauvcomm <- comm
    groups <- eval.parent(LAAA$groups)

    if(is.null(colnames(comm))) colnames(comm) <- paste("species", 1:ncol(comm), sep="")
    if(is.null(rownames(comm))) rownames(comm) <- paste("community", 1:nrow(comm), sep="")
    if(any(colSums(comm) < tol)){ 
        comm <- comm[, colSums(comm) > tol, drop=FALSE]
    }
    if(any(rowSums(comm) < tol)){ 
        comm2 <- comm[rowSums(comm) > tol, , drop=FALSE]
        groups <- groups[rowSums(comm) > tol]
        comm <- comm2
    }
    comm[comm < 0] <- 0
    if(!inherits(groups, "factor") | length(levels(groups)) != length(unique(groups)))
        groups <- factor(groups)

    method <- dbobj$method

    if(length(dbobj$test$pvalue) ==1){
        species <- FALSE
        global <- TRUE
        signif <- FALSE
    }
    else{
        species <- TRUE
        if(length(dbobj$test$pvalue) > ncol(comm))
        global <- TRUE
        else 
        global <- FALSE
    }
    if(method == "BrayCurtis" && signif){
        if(global) pval <- dbobj$test$pvalue[-1]
        else pval <- dbobj$test$pvalue
        if(all(pval > salpha)) stop("None of the species tests was significant")
        commE <- comm[, pval <= salpha]
        attributes(commE)$internal.p.rowSums <- rowSums(comm)
        comm <- commE
    }
    if(method != "BrayCurtis" && signif){
        if(global) pval <- dbobj$test$pvalue[-1]
        else pval <- dbobj$test$pvalue
        if(all(pval > salpha)) stop("None of the species tests was significant")
        comm <- comm[, pval <= salpha]
    }
    FF <- levels(groups)
    combFF <- combn(FF, 2)
    if(!species) {
        FUNcomb <- function(i){
        commC <- sauvcomm[groups%in%combFF[,i], ]
        groupsC <- factor(groups[groups%in%combFF[,i]])
        return(dbMANOVAspecies(commC, groupsC, nrep = nrep, method = method, global = TRUE, species = FALSE, padjust = padjust, tol = tol))
    }
    simuglobal <- lapply(1:ncol(combFF), FUNcomb)
    names(simuglobal) <- apply(combFF, 2, function(x) paste(x, collapse=":")) 
        simu <- simuglobal
        class(simu) <- "dbMANOVAspecies_pairwise"
        return(simu)
    }
    if(!signif | !global) {
    FUNcomb <- function(i){
        commC <- comm[groups%in%combFF[,i], ]
        groupsC <- factor(groups[groups%in%combFF[,i]])
        return(dbMANOVAspecies(commC, groupsC, nrep = nrep, method = method, global = global, species = TRUE, padjust = padjust, tol = tol))
    }
    simu <- lapply(1:ncol(combFF), FUNcomb)
    names(simu) <- apply(combFF, 2, function(x) paste(x, collapse=":")) 
    }
    else {
    FUNcomb <- function(i){
        commC <- comm[groups%in%combFF[,i], ]
        groupsC <- factor(groups[groups%in%combFF[,i]])
        return(dbMANOVAspecies(commC, groupsC, nrep = nrep, method = method, global = FALSE, species = TRUE, padjust = padjust, tol = tol))
    }
    simuspecies <- lapply(1:ncol(combFF), FUNcomb)
    names(simuspecies) <- apply(combFF, 2, function(x) paste(x, collapse=":")) 
    FUNcomb <- function(i){
        commC <- sauvcomm[groups%in%combFF[,i], ]
        groupsC <- factor(groups[groups%in%combFF[,i]])
        return(dbMANOVAspecies(commC, groupsC, nrep = nrep, method = method, global = TRUE, species = FALSE, padjust = padjust, tol = tol))
    }
    simuglobal <- lapply(1:ncol(combFF), FUNcomb)
    names(simuglobal) <- apply(combFF, 2, function(x) paste(x, collapse=":")) 
    simu <- list(global_test = simuglobal, per_species_test = simuspecies)
    }
    attributes(simu)$species.names <- colnames(sauvcomm)
    class(simu) <- "dbMANOVAspecies_pairwise"  
    return(simu)
}
