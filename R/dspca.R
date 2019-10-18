dspca <- function(comm, S=NULL, tol=1e-8){

    df <- t(comm)
    if(is.null(S)){
        S <- matrix(diag(rep(1, nrow(df)))) # By default: minimum similarity between any two species
        colnames(S) <- rownames(S) <- rownames(df)
    }
    if(!inherits(df, "data.frame"))
        df <- as.data.frame(df)
    dfp <- t(t(df)/colSums(df))
    step1 <- S
    svd.step1 <- svd(step1)
    u <- svd.step1$u
    d <- svd.step1$d
    r1 <- sum(d > (d[1] * tol))
    dp <- d[1:r1]
    up <- u[, 1:r1]
    X <- up%*%diag(sqrt(dp))
    Y <- t(dfp)%*%X
    rownames(X) <- rownames(df)
    colnames(X) <- paste("SPC", 1:ncol(X), sep="")
    colnames(Y) <- paste("SPC", 1:ncol(Y), sep="")
    Scom <- Y%*%t(Y)
    Q <- diag(1/sqrt(diag(Scom)))
    Y1 <- Q%*%Y
    rownames(Y1) <- rownames(Y)
    colnames(Y1) <- paste("SPC", 1:ncol(Y1), sep="")
    Y1 <- Q%*%Y
    Scom <- Q%*%Scom%*%Q
    rownames(Scom) <- colnames(Scom) <- colnames(df)
    step2 <- t(Y1)%*%Y1
    svd.step2 <- svd(step2)
    d <- svd.step2$d
    r2 <- sum(d > (d[1] * tol))
    dp <- d[1:r2]
    u <- svd.step2$u
    up <- u[, 1:r2]
    step2.Y <- Y1%*%up
    rownames(step2.Y) <- colnames(df)
    step2.X <- X%*%up
    rownames(step2.X) <- rownames(df)
    colnames(step2.X) <- paste("CPC", 1:ncol(step2.X), sep="")
    colnames(step2.Y) <- paste("CPC", 1:ncol(step2.Y), sep="")
    res <- list()
    res$eig <- dp
    res$X <- step2.X
    res$Y <- step2.Y
    res$Scom <- Scom
    class(res) <- "dspca"
    return(res)

}
