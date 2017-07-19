distMS <-
function (comm, diag = FALSE, upper = FALSE) 
{
 
    if (!(inherits(comm, "data.frame") | inherits(comm, "matrix"))) 
        stop("comm is not a data.frame or a matrix")
    nlig <- nrow(comm)
    d <- matrix(0, nlig, nlig)
    d.names <- row.names(comm)
    fun1 <- function(x) {
        sum(abs(comm[x[1], ] - comm[x[2], ]))/sum(apply(cbind.data.frame(comm[x[1], ], comm[x[2], ]), 1, max))
    }
    comm <- as.matrix(comm)
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])

    d <- unlist(apply(index, 1, fun1))

    attr(d, "Size") <- nlig
    attr(d, "Labels") <- d.names
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- "Marczewski-Steinhaus"
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
