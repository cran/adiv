DJac <- function (comm, dis, diag = FALSE, upper = FALSE) 
{
    if (!(inherits(comm, "data.frame") | inherits(comm, 
        "matrix"))) 
        stop("comm is not a data.frame or a matrix")
    if (!inherits(dis, "dist")) 
        stop("dis is not an object of class 'dist'")
    df <- as.matrix(comm)
    dis <- as.matrix(dis)
    if (nrow(dis) != ncol(df)) 
        stop("the size of dis is not equal to the number of columns in comm")
    if (is.null(colnames(df)) | is.null(rownames(dis))) 
        warning("Species (attributes) in dis was been assumed to be in the same order as species (columns) in comm")
    else {
        if (!all(c(colnames(df) %in% rownames(dis), rownames(dis) %in% 
            colnames(df)))) 
            warning("Species (attributes) in dis was been assumed to be in the same order as species (columns) in comm")
        else df <- df[, rownames(dis)]
    }
    if (!is.numeric(df)) 
        stop("df must contain  numeric values")
    if (any(df < 0)) 
        stop("non negative value expected in df")
    nlig <- nrow(df)
    nsp <- ncol(df)
    d.names <- row.names(df)
    if (is.null(d.names)) 
        d.names <- 1:nlig
    df <- as.matrix(1 * (df > 0))
    a <- df %*% t(df)
    b <- df %*% (1 - t(df))
    c <- (1 - df) %*% t(df)
    B <- C <- matrix(0, nlig, nlig)
    for (i in 2:nlig) {
        for (j in 1:(i - 1)) {
            compb <- (1:nsp)[df[i, ] > 1e-08 & df[j, ] < 1e-08]
            compY <- (1:nsp)[df[j, ] > 1e-08]
            bY <- dis[compb, compY]
            if (is.null(compb)) 
                B[i, j] <- 0
            else {
                if (length(compb) == 1) 
                  B[i, j] <- min(bY)
                else {
                  if (length(compY)==1)   B[i, j] <- sum(bY)
                  else B[i, j] <- sum(apply(bY, 1, min))
                }
            }
            compc <- (1:nsp)[df[j, ] > 1e-08 & df[i, ] < 1e-08]
            compX <- (1:nsp)[df[i, ] > 1e-08]
            cX <- dis[compc, compX]
            if (is.null(compc)) 
                C[i, j] <- 0
            else {
                if (length(compc) == 1) 
                  C[i, j] <- min(cX)
                else {
                  if (length(compX)==1)   C[i, j] <- sum(cX)
                  else    C[i, j] <- sum(apply(cX, 1, min))
                } 
            }
        }
    }
    Bf <- B + t(C)
    Cf <- C + t(B)
    B <- Bf
    C <- Cf
    A <- a + (b - B) + (c - C)
    j <- (B + C)/(a + b + c)
    jrich <- abs(B - C)/(a + b + c)
    jrepl <- j - jrich
    J <- as.dist(j)
    JRich <- as.dist(jrich)
    JRepl <- as.dist(jrepl)
    attr(J, "Size") <- attr(JRich, "Size") <- attr(JRepl, 
        "Size") <- nlig
    attr(J, "Labels") <- attr(JRich, "Labels") <- attr(JRepl, 
        "Labels") <- d.names
    attr(J, "Diag") <- attr(JRich, "Diag") <- attr(JRepl, 
        "Diag") <- diag
    attr(J, "Upper") <- attr(JRich, "Upper") <- attr(JRepl, 
        "Upper") <- upper
    attr(J, "method") <- "Jaccard"
    attr(JRich, "method") <- "Component of richness difference"
    attr(JRepl, "method") <- "Component of replacement or turnover"
    class(J) <- "dist"
    class(JRepl) <- "dist"
    class(JRich) <- "dist"
    return(list(J = J, JRepl = JRepl, JRich = JRich))
}

