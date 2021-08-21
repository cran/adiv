betaUniqueness <- function (comm, dis, Nind = 10000) 
{
    if (!(inherits(comm, "data.frame") | inherits(comm, "matrix"))) 
        stop("comm is not a data.frame or a matrix")
    if(nrow(comm) < 2) stop("comm must have at least two rows")
    D <- as.matrix(dis)
    if (any(D > 1)) {
        D <- D/max(D)
        warnings("All values in dis have been divided by their maximum (highest observed value)")
    }
    U <- dislptransport(comm, D, diag = FALSE, upper = FALSE, Nind = Nind)
    dis0 <- matrix(1, ncol(comm), ncol(comm))-diag(rep(1, ncol(comm))) 
    U0 <- dislptransport(comm, dis0, diag = FALSE, upper = FALSE, Nind = Nind)
    Ustar <- U0
    Ustar[Ustar<=0] <- 1
    res <- list()
    res$betaUniqueness <- as.matrix(U/Ustar)
    res$betaRedundancy <- 1 - as.matrix(U/Ustar)
    res$dissimilarityGap <- as.matrix(U0-U)
    res$DR <- as.matrix(U0)
    res$DKG <- as.matrix(U)
    return(res)
}
