treeUniqueness <- function(phyl, comm, index = c("richness", "GiniSimpson", "Shannon"), tol = 0.001)
{
    if (!inherits(comm, "data.frame")){
    if(!inherits(comm, "matrix"))
    stop("Object comm must be of class data.frame or matrix")
    }
    if (any(comm < 0))
        stop("Negative value in comm are not allowed")
    if (any(rowSums(comm) < 1e-16))
        stop("Empty plots in comm must be removed")
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
    if (is.ultrametric(phyl.phylo, tol = tol)== FALSE)
        warning("the phylogenetic phyl is not ultrametric")
    index <- index[1]
    if(!index%in%c("Shannon", "GiniSimpson", "richness"))
        stop("index must be Shannon or GiniSimpson or richness")
    TEC <- tecAptree(phyl.phylo)
    groups <- TEC$list
    index_base <- speciesdiv(comm, method = index)[, 1]
    funmerge <- function(w)
    {
       merged <- t(apply(comm, 1, function(x) tapply(x, w, FUN=sum)))
       return(speciesdiv(merged, method = index)[, 1])
    } 
    index_phyls <- lapply(groups, funmerge)
    all_branch_index <- t(t(cbind.data.frame(index_phyls)) * TEC$plength)
    mean_index_phyl <- rowSums(all_branch_index) / max(TEC$h)
    phylo_uniq <- mean_index_phyl / index_base
    phylo_redundancy <- 1 - phylo_uniq
    outputs <- data.frame(index_base,mean_index_phyl,phylo_uniq,phylo_redundancy)
    colnames(outputs) <- c("Dk", "Dp", "Phylogenetic Uniqueness", "Phylogenetic Redundancy")
    rownames(outputs) <- rownames(comm)
    return(outputs)
}
