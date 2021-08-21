dislptransport <-
function (comm, dis, diag = FALSE, upper = FALSE, Nind = 10000) 
{
 
    if (!(inherits(comm, "data.frame") | inherits(comm, "matrix"))) 
        stop("comm is not a data.frame or a matrix")
    if(nrow(comm) < 2) stop("comm must have at least two rows")
    D <- as.matrix(dis)
    nlig <- nrow(comm)
    ncol <- ncol(comm)
    d <- matrix(0, nlig, nlig)
    d.names <- row.names(comm)
    df <- as.data.frame(t(comm))
    dfp <- t(t(df)/colSums(df))
    funlp <- function(x) {
        row_rhs <- round(dfp[,x[1]] * Nind)
        if(sum(row_rhs) < Nind) {
            nmiss <- Nind - sum(row_rhs)
            e <- sample((1:nrow(dfp))[row_rhs < dfp[,x[1]] * Nind], nmiss)
            for(i in e)    row_rhs[i] <- row_rhs[i] + 1
        }
        if(sum(row_rhs) > Nind) {
            ntoo <- abs(Nind - sum(row_rhs))
            e <- sample((1:nrow(dfp))[row_rhs > dfp[,x[1]] * Nind], ntoo)
            for(i in e)    row_rhs[i] <- row_rhs[i] - 1
        }
        col_rhs <- round(dfp[,x[2]] * Nind)
        if(sum(col_rhs) < Nind) {
            nmiss <- Nind - sum(col_rhs)
            e <- sample((1:nrow(dfp))[col_rhs < dfp[,x[2]] * Nind], nmiss)
            for(i in e)    col_rhs[i] <- col_rhs[i] + 1
        }
        if(sum(col_rhs) > Nind) {
            ntoo <- abs(Nind - sum(col_rhs))
            e <- sample((1:nrow(dfp))[col_rhs > dfp[,x[2]] * Nind], ntoo)
            for(i in e)    col_rhs[i] <- col_rhs[i] - 1
        }
        row_signs <- rep("<=", ncol)
        col_signs <- rep(">=", ncol)
        res <- sum(lp.transport(D, "min", row_signs, row_rhs, col_signs, col_rhs)$solution*D/Nind)
       return(res)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    d <- unlist(apply(index, 1, funlp))
    attr(d, "Size") <- nlig
    attr(d, "Labels") <- d.names
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- "GregoriusGilletZiehe"
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
