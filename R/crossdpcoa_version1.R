crossdpcoa_version1 <-
function(df, facA, facB, dis = NULL, scannf
= TRUE, nf = 2, w = c("classic", "independence"), tol = 1e-07) {
#######################################################
#                           Check                     #
#######################################################

    if (!inherits(df, "data.frame") & !inherits(df, "matrix")) 
        stop("df must be a data frame or a matrix")
    df <- t(df)
    if (!is.factor(facA))
        stop("facA is not a factor")
    if (!is.factor(facB))
        stop("facB is not a factor")
    r.n <- row.names(df)
    if (any(df < -tol))
        stop("negative entries in df")
    if (any(df < 0))
       df[df < 0] <- 0


################################################################
#                      Preparing data                          #
################################################################

    nsp <- length(r.n)
    nA <- length(levels(facA))
    nB <- length(levels(facB))
    nC <- ncol(df)


    if(is.null(dis)){
        dis <- (matrix(1, nsp, nsp) - diag(rep(1, nsp))) *
            sqrt(2)
        rownames(dis) <- colnames(dis) <- r.n
        dis <- as.dist(dis)
    }

    # Vectors of weights

    ntot <- sum(df)
    wij <- colSums(df)/sum(df)
    lambdai <- tapply(wij, facA, sum)
    muj <- tapply(wij, facB, sum)
    names(lambdai) <- levels(facA)
    longlambdai <- lambdai[as.character(facA)]
    names(muj) <- levels(facB)
    longmuj <- muj[as.character(facB)]
    PbC <- sweep(df, 2, colSums(df), "/")

    if(w[1] == "independence"){
        wij <- longlambdai*longmuj
        epsilonk <- apply(PbC, 1, function(x) sum(wij * x))
    } else if(w[1] == "classic"){
        epsilonk <- rowSums(df)/ntot
    } else if(is.numeric(w)){
        wij <- w / sum(w)
        epsilonk <- apply(PbC, 1, function(x) sum(wij * x))
        lambdai <- tapply(wij, facA, sum)
        muj <- tapply(wij, facB, sum)
        names(lambdai) <- levels(facA)
        longlambdai <- lambdai[as.character(facA)]
        names(muj) <- levels(facB)
        longmuj <- muj[as.character(facB)]
    }
    else stop("Incorrect definition of parameter w")
    PbA <- apply(PbC, 1, function(x) tapply(x*wij/longlambdai, facA, sum))
    PbA <- as.data.frame(t(PbA))
    PbB <- apply(PbC, 1, function(x) tapply(x*wij/longmuj, facB, sum))
    PbB <- as.data.frame(t(PbB))

#################################################################
#                          Computations                         #
#################################################################

    pcosp <- dudi.pco(dis, row.w = epsilonk, scannf = FALSE, full = TRUE, tol = tol)
    X <- as.matrix(pcosp$li)
    Bcat <- diag(epsilonk)
    YA <- t(PbA)%*%X
    YB <- t(PbB)%*%X
    YC <- t(PbC)%*%X
    BA <- diag(lambdai)
    BB <- diag(muj)
    BC <- diag(wij)
    BX <- diag(epsilonk)

    # We are in the space of APQE
    # Calculation of the components of diversity:

    DIVA <- sum(colSums(BA%*%((YA)^2)))

    DIVB <- sum(colSums(BB%*%((YB)^2)))

    DIVXtotale <- sum(colSums(BX%*%((X)^2)))

    C <- as.matrix(df)
    sumcolC <- diag(1/colSums(C))
    Cp <- C%*%sumcolC
    meanCp <- t(Cp)%*%X
    sumCCp <- t(Cp)%*%X^2
    varCCp <- rowSums(sumCCp-meanCp^2)
    DIVintraC <- sum(wij*varCCp)

    DIVAB <- DIVXtotale - (DIVintraC + DIVA + DIVB)
    divcomponents <- c(DIVintraC, DIVA, DIVB, DIVAB, DIVXtotale)

# Recentring process: the new positions of the levels of factor B are all in the center of the space.

    mB <- model.matrix(~-1+facB)
    projB <- mB%*%solve(t(mB)%*%diag(wij)%*%mB)%*%t(mB)%*%diag(wij)
    projBo <- diag(rep(1, nC)) - projB

    XaB <- X
    YAaB <- YA
    YCaB <- projBo%*%YC

# In this space, the principal axes of factor A are defined and all points are projected on these axes

    svdYAaB <- svd(t(YAaB)%*%BA%*%YAaB)
    XoAaB <- XaB%*%svdYAaB$u
    YAoAaB <- YAaB%*%svdYAaB$u
    YCoAaB <- YCaB%*%svdYAaB$u
    colnames(XoAaB) <- colnames(YAoAaB) <- colnames(YCoAaB) <- paste("Axis", 1:ncol(YCoAaB), sep="")

    eig <- svdYAaB$d[svdYAaB$d>tol]
    rank <- length(svdYAaB$d[svdYAaB$d>tol])

    if (scannf) {
            barplot(eig[1:rank])
            cat("Select the number of axes: ")
            nf <- as.integer(readLines(n = 1))
    }
    if (nf < 1)
        nf <- 2
    if(nf > rank) nf <- rank
    names(divcomponents) <- c("SSW", "SSA", "SSB", "SSAB", "SST")
    res <- list()

    rownames(YCoAaB) <- colnames(df)
    res$l1 <- as.data.frame(XoAaB[, 1:nf])
    res$l2 <- as.data.frame(YAoAaB[, 1:nf])
    res$l3 <- as.data.frame(YCoAaB[, 1:nf])
    res$eig <- eig
    res$lwX <- epsilonk
    res$lwA <- lambdai
    res$lwB <- muj
    res$lwC <- wij

    res$div <- divcomponents

    res$call <- match.call()
    return(res)

}
