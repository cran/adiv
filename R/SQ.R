SQ <-
function(comm, Sigma = NULL, type=c("similarity", "dissimilarity")){

    type <- type[1]
    if(!type%in%c("similarity", "dissimilarity"))
        stop("Incorrect definition of parameter type")
    if (!is.null(Sigma)){
    if (type=="similarity")
    return(1-as.matrix(disc(as.data.frame(t(comm)), as.dist(sqrt(2*(1-Sigma)))))^2/2)
    else
     return(disc(as.data.frame(t(comm)), as.dist(sqrt(2*(1-Sigma))))^2/2)
    }
    else{
    if (type=="similarity")
    return(1-as.matrix(disc(as.data.frame(t(comm))))^2/2)
    else
    return(disc(as.data.frame(t(comm)))^2/2)
    }
}
