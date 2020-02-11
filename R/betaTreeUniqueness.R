betaTreeUniqueness <- function(mtree, comm, height = NULL, tol = 0.001){

    BrayCurtis <- function (COMM, TOL) 
    {
    COMM <- data.frame(COMM)
    if (!inherits(COMM, "data.frame")) 
        stop("COMM is not a data.frame")
    nrow <- nrow(COMM)
    d <- matrix(0, nrow, nrow)
    dnames <- row.names(COMM)
    FUN <- function(x) {
        sum(abs(COMM[x[1], ] - COMM[x[2], ]))/ sum(COMM[x[1], ] + COMM[x[2], ])
    }
    COMM <- as.matrix(COMM)
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, FUN))

    attr(d, "Size") <- nrow
    attr(d, "Labels") <- dnames
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    class(d) <- "dist"
    return(d)
    }

    Dp <- DP(mtree = mtree, comm = comm, height = height, 
        diag = FALSE, upper = FALSE, tol = tol)

    Ds <- BrayCurtis(comm)

    Up <- as.matrix(Dp/Ds)

    return(Up)

}


