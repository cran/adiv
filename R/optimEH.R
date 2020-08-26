optimEH <-
function(phyl, nbofsp, tol = 1e-8, give.list = TRUE)
{
    arg.phyl <- .checkphyloarg(phyl)
    phyl <- arg.phyl$phyl
    phyl.phylo <- arg.phyl$phyl.phylo
    rm(arg.phyl)

    if(is.binary.phylo(phyl.phylo))
        phy.h <- as.hclust(phyl.phylo) ## also test for ultrametricity 
    else{
        if(!is.ultrametric(phyl.phylo)) stop("the tree is not ultrametric")
        phy.h <- hclust(as.dist(cophenetic.phylo(phyl.phylo)/2), "average")
    }
    phyl.D <- cophenetic.phylo(phyl.phylo)/2
    nbesp <- nTips(phyl)
  
    if (length(nbofsp) != 1) stop("unconvenient nbofsp")
    nbofsp <- round(nbofsp)
    if (!((1 < nbofsp) & (nbofsp <= nbesp))) stop("nbofsp must be between 2 and ", nbesp)
    sp.names <- phy.h$labels
    if (nbofsp == nbesp) {
        res1 <- EH(phyl)
        sauv.names <- sp.names
    }
    else {
        phyl.D <- as.matrix(phyl.D)
        num.Orig <- as.vector(solve(phyl.D, rep(1, nbesp)))
        denum.Orig <- as.vector(t(rep(1, nbesp))%*%num.Orig)
        Orig <- as.data.frame(num.Orig/denum.Orig)
        rownames(Orig) <- rownames(phyl.D)
        car1 <- split(Orig, cutree(phy.h, nbofsp))
        name1 <- lapply(car1,function(x) rownames(x)[abs(x - max(x)) < tol])
        sauv.names <- lapply(name1, paste, collapse = " OR ")
        comp <- as.character(as.vector(lapply(name1, function(x) x[1])))
       res1 <- EH(phyl, select = comp)
    }
    if (give.list == TRUE)
        return(list(value = res1, selected.sp = cbind.data.frame(names = unlist(sauv.names))))
    else
        return(res1)
}
