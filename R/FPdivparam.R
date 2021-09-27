FPdivparam <- function (comm, disORtree, method = c("KY", "KstarI"), palpha = 2, equivalent = FALSE, option = c("asymmetric", "symmetric"), dmax = NULL, tol = 1e-8) 
{

    option <- option[1]
    method <- method[1]
    if (!method %in% c("KY", "KstarI")) 
        stop("unconvenient method")
    if (!inherits(comm, "data.frame") & !inherits(comm, 
        "matrix")) 
        stop("comm must be a data frame or a matrix")
    comm <- as.matrix(comm)
    if (any(comm < 0)) 
        stop("Negative value in comm")
    nsp <- ncol(comm)
    if (!is.null(disORtree)) {
        if(inherits(disORtree, "phylo") | inherits(disORtree, "phylo4") | inherits(disORtree, "hclust")){
        arg.phyl <- .checkphyloarg(disORtree)
        phyl.phylo <- arg.phyl$phyl.phylo
        tre4 <- arg.phyl$phyl
        if (!hasEdgeLength(tre4))
            phyl.phylo <- compute.brlen(phyl.phylo, 1)
        if(is.ultrametric(phyl.phylo) | option == c("symmetric"))
            dis <- cophenetic.phylo(phyl.phylo)/2
        else
            dis <- matrix(rep(diag(vcv.phylo(phyl.phylo))), nsp, nsp)-vcv.phylo(phyl.phylo)
        }
        else if(inherits(disORtree, "matrix") | inherits(disORtree, "dist") )
         dis <- as.matrix(disORtree)
    else stop("unconvenient definition of disORtree")
        if(any(!colnames(comm)%in%colnames(dis))) 
        stop("Some species names in comm are missing in disORtree")
        dis <- dis[colnames(comm), colnames(comm)]
        dis <- as.dist(dis)
    }
    if (is.null(disORtree)) {
        dis <- as.dist( matrix(1, ncol(comm), ncol(comm)) 
            - diag(rep(1, ncol(comm))) )
    }

if(method == "KY") {
    if(max(dis) > 1) {
        if(!is.null(dmax) && (dmax-max(dis)) < tol)  dis <- dis/dmax
        else dis <- dis/max(dis)
    }
    sim <- 1 - as.matrix(dis)
    FREQ <- sweep(comm, 1, rowSums(comm), "/")
    FUN <- function(aaa) {
       divv <- rep(0, nrow(comm))
    for (i in 1:nrow(comm)) {
        if (sum(comm[i, ]) < 1e-16) 
            divv[i] <- 0
        else{
            FFF <- FREQ[i, FREQ[i, ]> tol]
            simFFF <- sim[FREQ[i, ]> tol, FREQ[i, ]> tol]
            if(abs(aaa-1) > tol){
            if(equivalent)
            divv[i] <- (t(FFF) %*% t(((FFF%*% simFFF)^(aaa-1)) ))^(1/(1-aaa))
            else
            divv[i] <- t(FFF) %*% t((1-(FFF%*% simFFF)^(aaa-1)) )/(aaa-1)
            }
            else{
                divv[i] <- -t(FFF) %*% t(log(FFF%*% simFFF))
                if(equivalent) divv[i] <- exp(divv[i])
            }
            }
    }
        return(divv)
    }
    if(length(palpha)==1){
       div <- cbind.data.frame(FUN(palpha))
       colnames(div) <- "K"
       rownames(div) <- rownames(comm)
       class(div) <- c("FPdivparam", "data.frame")
    }
    else{
    div <- sapply(palpha, FUN)
    colnames(div) <- paste("K", 1:length(palpha), sep="_")
    rownames(div) <- rownames(comm)
    div <- list(palpha = palpha, div = div)
    class(div) <- c("FPdivparam", "list")
    }
    return(div)
 }
 else{
    commgardees <- (1:nrow(comm))[rowSums(comm) >= tol]
    commdiv <- comm[commgardees, ]
    FREQ <- sweep(commdiv, 1, rowSums(comm), "/")
    PA <- commdiv
    PA[PA>0] <- 1

    FUN <- function(aaa) {
       ORIval <- distinctAb(comm = commdiv, disORtree = dis, method = "KstarI", palpha = aaa)[[1]]
       ORIval[is.na(ORIval)] <- 0
       divv <- rowSums(ORIval)
       if(equivalent){
       if(is.null(dmax))  dmax <- max(dis)
       if((dmax-max(dis)) >= tol)   dmax <- max(dis) 
       if(abs(aaa-1) > tol) divv <- (( 1-(aaa-1)*divv/dmax ))^(1/(1-aaa))
       else divv <- exp(divv/dmax)
       }
       return(divv)
    }
    if(length(palpha)==1){
       div <- FUN(palpha)
       names(div) <- rownames(commdiv)
       DIV <- rep(0, nrow(comm))
       names(DIV) <- rownames(comm)
       DIV[names(div)] <- div
       class(div) <- c("FPdivparam", "vector")

    }
    else{
        div <- sapply(palpha, FUN)
        rownames(div) <- rownames(commdiv)
        DIV <- matrix(0, nrow(comm), length(palpha))
        rownames(DIV) <- rownames(comm)
        colnames(DIV) <- paste("KstarI", 1:length(palpha), sep="_")
        DIV[rownames(div), ] <- div
        div <- DIV
        div <- list(palpha = palpha, div = DIV)
        class(div) <- c("FPdivparam", "list")
    }
    return(div)
    }
}
