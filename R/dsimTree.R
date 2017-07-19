dsimTree <-
function(phyl, method = c(1, 2, 3, 4, 5),
rootedge = NULL, type=c("similarity", "dissimilarity")){

    type <- type[1]
    if(!type%in% c("dissimilarity", "similarity")) stop("type must be either dissimilarity or similarity")
    tree <- .checkphyloarg(phyl)$phyl.phylo
    method <- method[1]
    a <- ape::vcv(tree)
    if(!is.null(rootedge))
    a <- a + rootedge
    ab <- diag(a)%*%t(rep(1, ncol(ape::vcv(tree))))
    ac <- rep(1, ncol(ape::vcv(tree)))%*%t(diag(a))
    b <- ab-a
    c <- ac-a
    if(type=="similarity"){
    # taxonomic similarities among species
    if (method==1)
    return(a/(a+2*b+2*c))
    if (method==2)
    return(a/(a+b+c))
    if (method==3)
    return(2*a/(2*a+b+c))
    if (method==4)
    return(a/sqrt((a+b)*(a+c)))
    if (method==5)
    return(4*a/(4*a+b+c))
    }
    else{
    # taxonomic disimilarities among species
    if (method==1)
    return(as.dist(1-a/(a+2*b+2*c)))
    if (method==2)
    return(as.dist(1-a/(a+b+c)))
    if (method==3)
    return(as.dist(1-2*a/(2*a+b+c)))
    if (method==4)
    return(as.dist(1-a/sqrt((a+b)*(a+c))))
    if (method==5)
    return(as.dist(1-4*a/(4*a+b+c)))

    }

}
