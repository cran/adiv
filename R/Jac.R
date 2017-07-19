Jac <-
function (comm, diag = FALSE, upper = FALSE) 
{
    if (!(inherits(comm, "data.frame") | inherits(comm, "matrix"))) 
        stop("comm is not a data.frame or a matrix")
    df <- as.matrix(comm)
    if (!is.numeric(df)) 
        stop("df must contain  numeric values")
    if (any(df < 0)) 
        stop("non negative value expected in df")
    nlig <- nrow(df)
    d.names <- row.names(df)
    if (is.null(d.names)) 
        d.names <- 1:nlig
    df <- as.matrix(1 * (df > 0))
    a <- df %*% t(df)
    b <- df %*% (1 - t(df))
    c <- (1 - df) %*% t(df)
    j <- (b+c)/(a + b + c)
    jrich <- abs(b-c)/(a+b+c)
    jrepl <- j-jrich
    J <- as.dist(j)
    JRich <- as.dist(jrich)
    JRepl <- as.dist(jrepl)
    attr(J, "Size") <- attr(JRich, "Size") <- attr(JRepl, "Size") <- nlig
    attr(J, "Labels") <- attr(JRich, "Labels") <- attr(JRepl, "Labels") <- d.names
    attr(J, "Diag") <- attr(JRich, "Diag") <- attr(JRepl, "Diag") <- diag
    attr(J, "Upper") <- attr(JRich, "Upper") <- attr(JRepl, "Upper") <-  upper
    attr(J, "method") <- "Jaccard"
    attr(JRich, "method") <- "Component of richness difference"
    attr(JRepl, "method") <- "Component of replacement or turnover"
    class(J) <- "dist"
    class(JRepl) <- "dist"
    class(JRich) <- "dist"
    return(list(J=J, JRepl = JRepl, JRich = JRich))
}
