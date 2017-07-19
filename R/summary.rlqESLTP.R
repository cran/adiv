summary.rlqESLTP <-
function (object, ...) 
{
    if (!inherits(object, "rlqESLTP")) 
        stop("to be used with 'rlqESLTP' object")
    appel <- as.list(object$call)
    dudiL <- object$dudiL
    dudiR <- object$dudiR
    dudiQ <- object$dudiQ
    norm.w <- function(object, w) {
        f2 <- function(v) sqrt(sum(v * v * w)/sum(w))
        norm <- apply(object, 2, f2)
        return(norm)
    }
    util <- function(n) {
        x <- "1"
        for (i in 2:n) x[i] <- paste(x[i - 1], i, sep = "")
        return(x)
    }
    eig <- object$eig[1:object$nf]
    covar <- sqrt(eig)
    sdR <- norm.w(object$lR, dudiR$lw)
    sdQ <- norm.w(object$lQ, dudiQ$lw)
    corr <- covar/sdR/sdQ
    U <- cbind.data.frame(eig, covar, sdR, sdQ, corr)
    row.names(U) <- as.character(1:object$nf)
    cat("\nEigenvalues decomposition:\n")
    print(U)
    cat("\nInertia & coinertia R:\n")
    inertia <- cumsum(sdR^2)
    max <- cumsum(dudiR$eig[1:object$nf])
    ratio <- inertia/max
    U <- cbind.data.frame(inertia, max, ratio)
    row.names(U) <- util(object$nf)
    print(U)
    cat("\nInertia & coinertia Q:\n")
    inertia <- cumsum(sdQ^2)
    max <- cumsum(dudiQ$eig[1:object$nf])
    ratio <- inertia/max
    U <- cbind.data.frame(inertia, max, ratio)
    row.names(U) <- util(object$nf)
    print(U)
    cat("\nCorrelation L:\n")
    max <- sqrt(dudiL$eig[1:object$nf])
    ratio <- corr/max
    U <- cbind.data.frame(corr, max, ratio)
    row.names(U) <- 1:object$nf
    print(U)

}
