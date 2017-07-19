CFprop <-
function (df, wA=rep(1,ncol(df)))
{
    if (!inherits(df, "data.frame") & !inherits(df, "matrix")) 
        stop("df must be a data frame or a matrix")
    if (any(df < 0))
        stop("non-negative values expected in df")
    dfs <- apply(df, 1, sum)
    if (any(dfs == 0))
        stop("row with all zero value")
    if (any(wA < 0))
        stop("non-negative values expected in wA")
    df <- df/dfs
    nlig <- nrow(df)
    Cmat <- matrix(0, nlig, nlig)
    C.names <- row.names(df)
    df <- as.matrix(df)
    fun1 <- function(x) {
        p <- df[x[1], ]
        q <- df[x[2], ]
        w <- sum(wA * p * q)/sum(wA * p * p)/sum(wA * q * q)
        return(w)
    }
    index <- cbind(col(Cmat)[col(Cmat) < row(Cmat)], row(Cmat)[col(Cmat) < row(Cmat)])
    C <- unlist(apply(index, 1, fun1))
    for(i in 1:nrow(index))
    Cmat[index[i,1], index[i,2]] <- C[i]
    Cmat <- Cmat + t(Cmat)
    index <- cbind(1:nlig, 1:nlig)
    Cmat <- Cmat + diag(unlist(apply(index, 1, fun1)))
    rownames(Cmat)<-colnames(Cmat)<-C.names
    return(Cmat)
}
