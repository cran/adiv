dsimTaxo <-
function(tax, method = c(1, 2, 3, 4, 5), type=c("similarity", "dissimilarity")){

    type <- type[1]
    if(!type%in% c("dissimilarity", "similarity")) stop("type must be either dissimilarity or similarity") 
    if(!inherits(tax, "taxo")) stop("Object tax must be of class taxo")
    dtax <- dist.taxo(tax)^2
    tree <- hclust(dtax, method="average")
    tree <- as.phylo(tree)
    method <- method[1]
    a <- vcv.phylo(tree)
    ab <- diag(a)%*%t(rep(1, ncol(a)))
    ac <- rep(1, ncol(a))%*%t(diag(a))
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
