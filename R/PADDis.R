PADDis <-
function (comm, dis, method = NULL, diag = FALSE, upper = FALSE) 
{
    METHODS <- c("components", "Jaccard", "Sorensen", "SockalSneath", "Ochiai", "Simpson")
    if (!(inherits(comm, "data.frame") | inherits(comm, "matrix"))) 
        stop("comm is not a data.frame or a matrix")
    if (!inherits(dis, "dist"))
        stop("dis is not an object of class 'dist'")
    df <- as.matrix(comm)
    dis <- as.matrix(dis)
    if(nrow(dis)!=ncol(df))
        stop("the size of dis is not equal to the number of columns in comm")
    if(is.null(colnames(df)) | is.null(rownames(dis)))
	   warning("Species (attributes) in dis was been assumed to be in the same order as species (columns) in comm")
    else {if (!all(c(colnames(df)%in%rownames(dis), rownames(dis)%in%colnames(df))))
	   warning("Species (attributes) in dis was been assumed to be in the same order as species (columns) in comm")
    else
	    df <- df[, rownames(dis)]}
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
    if (is.null(method)) {
        cat("0 = list of all components (a, b, c, A, B, C)\n")
        cat("1 = Jaccard index (1901)\n")
        cat("d1 = (B+C)/(a+b+c)\n")
        cat("2 = Czekanowski (1913) or Sorensen (1948)\n")
        cat("d2 = (B+C)/(2*a+b+c)\n")
        cat("3 = Sockal & Sneath (1963)\n")
        cat("d3 = 2(B+C)/(a+2(b+c))\n")
        cat("4 = Ochiai (1957)\n")
        cat("d4 = [sqrt((A+B)(A+C))-A]/sqrt((a+b)(a+c))\n")
        cat("5 = Simpson (1943) \n")
        cat("d5 = min(B,C)/(a+min(b,c))\n")
        cat("6 = Kulczynski (1943) \n")
        cat("d6 = 0.5*(B/(a+b)+C/(a+c))\n")
        cat("Select an integer (0-6): ")
        method <- as.integer(readLines(n = 1))
    }
    a <- df %*% t(df)
    b <- df %*% (1 - t(df))
    c <- (1 - df) %*% t(df)
    B <- C <- matrix(0, nlig, nlig)
    for (i in 2:nlig){
	for(j in 1:(i-1)){
           compb <- (1:nsp)[df[i,]>1e-8 & df[j, ]<1e-8]
           compY <- (1:nsp)[df[j, ]>1e-8]
           bY <- dis[compb, compY]
           if(is.null(compb)) B[i, j] <- 0
           else {if(length(compb)==1) B[i, j] <- min(bY)
           else if(length(compY)==1) B[i, j] <- sum(bY)
           else	B[i, j] <- sum(apply(bY, 1, min)) }
           compc <- (1:nsp)[df[j,]>1e-8 & df[i, ]<1e-8]
           compX <- (1:nsp)[df[i, ]>1e-8]
           cX <- dis[compc, compX]
           if(is.null(compc)) C[i, j] <- 0
           else {if(length(compc)==1) C[i, j] <- min(cX)
           else if(length(compY)==1) C[i, j] <- sum(cX)
    		else C[i, j] <- sum(apply(cX, 1, min)) }
    }
    }
    Bf <- B + t(C)
    Cf <- C + t(B)
    B <- Bf
    C <- Cf
    colnames(B) <- colnames(C) <- rownames(B) <- rownames(C) <- rownames(df)
    A <- a + (b - B) + (c - C)
    if (method == 0) {
        return(list(a = a, b = b, c = c, A = A, B = B, C = C))
    }
    if (method == 1) {
        d <- (B+C)/(a + b + c)
    }
    else if (method == 2) {
        d <- (B+C)/(2 * a + b + c)
    }
    else if (method == 3) {
        d <- (2*(B+C))/(a + 2 * (b + c))
    }
    else if (method == 4) {
        d <- (sqrt(A+B)*sqrt(A+C)-A)/sqrt((a + b) * (a + c))
    }
    else if (method == 5) {
        minBC <- (B<=C)*B+(B>C)*C
        minbc <- (b<=c)*b+(b>c)*c
        d <- minBC / (a + minbc)
    }
    else if (method == 6) {
        d <- 0.5*(B/(a+b) + C/(a+c))
    }
    else stop("Non convenient method")
    d <- as.dist(d)
    attr(d, "Size") <- nlig
    attr(d, "Labels") <- d.names
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method+1]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
