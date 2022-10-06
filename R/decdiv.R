decdiv <- function(phyl, comm, dis = NULL, tol = 1e-08, option = 1:5, formula = c("QE", "EDI")){

    df <- comm
    if(!formula[1]%in%c("QE", "EDI")) stop("formula must be 'QE' or 'EDI'")
    if(formula[1]=="QE" & !is.null(dis)) dis <- sqrt(2*dis)
    if(!option[1]%in%(1:5)) stop("option must be 1, 2, 3, 4 or 5")
    if(is.vector(df)){
        if(is.null(names(df))) stop("comm must have names")
        vf <- df
        df <- cbind.data.frame(t(df))
        rownames(df) <- "comm1"
        colnames(df) <- names(vf)
    }
    if(!is.data.frame(df)) stop("df should be a data frame")
    if(any(rowSums(df)<tol)) stop("null row sums in comm: empty communities")
    if(any(df < -tol)) stop("negative values in df")
    df[df < tol] <- 0

    df <- as.data.frame(t(df))
    
    nsp <- nrow(df)

    if(!inherits(phyl, "phylo4"))
        tre4 <- .checkphyloarg(phyl)$phyl
    else tre4 <- phyl
    nsp <- length(tipLabels(tre4))

    if(!isRooted(tre4)){
        treape <- as(tre4, "phylo")
        treape$root.edge <- 0
        tre4 <- as(treape, "phylo4")
    }
    if(!hasNodeLabels(tre4)) nodeLabels(tre4) <- names(nodeLabels(tre4))
    else{
        e <- nodeLabels(tre4)
        e[is.na(e)] <- names(e[is.na(e)])
        nodeLabels(tre4) <- e
    }
    lTips <- descendants(tre4, getNode(tre4, type="internal"), type="tips")
    names(lTips) <- as.vector(nodeLabels(tre4))
    if(any(!rownames(df)%in%tipLabels(tre4)))
        stop("some species labels in comm are not in phyl")
    df <- df[tipLabels(tre4), , drop=FALSE]

    if (is.null(dis)){
        dis <- as.dist((matrix(1, nsp, nsp) - diag(rep(1,
            nsp))) * sqrt(2))
        attributes(dis)$Labels <- tipLabels(tre4)
    }

    if(!inherits(dis, "dist")) stop("dis must be of class dist")
    if(is.null(attributes(dis)$Labels))
        stop("dis must have labels")
    if(any(!attributes(dis)$Labels%in%rownames(df)))
        stop("some species labels in dis are not in comm")

    dis <- as.dist(as.matrix(dis)[rownames(df), rownames(df)])
    if(option[1]%in%(2:3) & (max(dis)^2/2)>1){
        DDis <- dis^2/2
        DDis <- DDis/max(DDis)
        dis <- sqrt(2*DDis)
    }
    lDD <- descendants(tre4, getNode(tre4, type="internal"), type = "children")
    names(lDD) <- as.vector(nodeLabels(tre4))
    mB <- cbind.data.frame(lapply(lTips, function(x) match(tipLabels(tre4), names(x), nomatch=0)))
    mB[mB > 0] <- 1
    rownames(mB) <- tipLabels(tre4)
    lDD <- lDD[colnames(mB)]

    diagesp <- diag(rep(1, nsp))
    rownames(diagesp) <- colnames(diagesp) <- tipLabels(tre4)
    mB <- cbind.data.frame(mB, diagesp)


    decdivV <- function(v) {

        mBab <- mB*v
        wnodes <- colSums(mBab[, 1:length(lDD)])/sum(v)
        names(wnodes) <- colnames(mBab[, 1:length(lDD)])
        fundiv <- function(x){
            if(sum(x) > 0){
            x <- x/sum(x)
            return(t(x)%*%as.matrix(dis^2/2)%*%x)
            }
            else return(0)
        }
        divnodes <- sapply(mBab[, 1:length(lDD)], fundiv)
        names(divnodes) <- colnames(mBab[, 1:length(lDD)])
        divnodes2 <- c(divnodes, rep(0, nsp))
        names(divnodes2) <- c(names(divnodes), tipLabels(tre4))
        wnodes2 <- c(wnodes, rep(0, nsp))
        names(wnodes2) <- c(names(wnodes), tipLabels(tre4))
        fundivintra <- function(i){
            return(sum(divnodes2[names(lDD[[i]])]*wnodes2[names(lDD[[i]])]))
        }
        divintra <- sapply(1:length(lDD), fundivintra)

        funM <- function(v){
            return(length(names(v)[colSums(mBab[, 1:length(lDD)])[names(v)]>tol]))
        }
        M <- as.vector(unlist(lapply(lDD, funM)))
        w <- wnodes[names(lDD)]
        w[w < tol] <- 1
        divintra <- divintra/w
        if(option[1] == 1)
            divfinal <- wnodes * (divnodes-divintra)
        else if(option[1] == 2)
            divfinal <- (divnodes-divintra) / (1-divintra) * M / (M - 1)
 	   else if(option[1] == 3)
            divfinal <- (divnodes-divintra) / (1-divnodes) / (M - 1) 
 	   else if(option[1] == 4)
            divfinal <- (divnodes-divintra)
        else
            divfinal <- wnodes
        return(divfinal)

    }

    res <- apply(df, 2, decdivV)
    attributes(res)$phyl <- as(tre4, "phylo")
    if(!is.null(colnames(df)))
        colnames(res) <- colnames(df)
    class(res) <- c("decdiv", "matrix")
    return(res)
    
}

