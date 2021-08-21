dbMANOVAspecies <- function(comm, groups, nrep = 999, method = c("Euclidean", "Manhattan", "Canberra", "BrayCurtis"), global = TRUE, species = TRUE, padjust = "none", tol = 1e-8){

    global <- global[1]
    species <- species[1]
    if(!inherits(global, "logical")) stop("Incorrect definition of argument global")
    if(!inherits(species, "logical")) stop("Incorrect definition of argument species")
    if(!species) global <- TRUE
    if(!inherits(comm, "data.frame") && !inherits(comm, "matrix")) stop("comm must be a data frame or a matrix")
    if(is.null(colnames(comm))) colnames(comm) <- paste("species", 1:ncol(comm), sep="")
    if(is.null(rownames(comm))) rownames(comm) <- paste("community", 1:nrow(comm), sep="")
    if(any(comm < -tol)) stop("comm must have nonnegative values")
    if(any(colSums(comm) < tol)){ 
        comm <- comm[, colSums(comm) > tol, drop=FALSE]
        warning("species with zero abundance over all considered communities have been removed")
    }
    if(nrow(comm) != length(groups)) stop("The length of argument groups should be equal to the number of rows of argument comm")
    if(any(rowSums(comm) < tol)){ 
        comm2 <- comm[rowSums(comm) > tol, , drop=FALSE]
        groups <- groups[rowSums(comm) > tol]
        comm <- comm2
        warning("empty communities with zero species abundance have been removed")
    }
    comm[comm < 0] <- 0
    if(!inherits(groups, "factor") && !inherits(groups, "character")) stop("Incorrect definition for argument groups")
    if(!inherits(groups, "factor") | length(levels(groups)) != length(unique(groups)))
        groups <- factor(groups)
    method <- method[1]
    if(!method%in%c("Euclidean","Manhattan","Canberra","BrayCurtis")) stop("Incorrect definition of argument method")
    Nplots <- nrow(comm)


if(!species){

distabund <- function (abund, meth) 
{

    d <- matrix(0, nrow(abund), nrow(abund))
    funEuclidean <- function(x) {
        sum((abund[x[1], ] - abund[x[2], ])^2)
    }
    funManhattan <- function(x) {
        sum(abs(abund[x[1], ] - abund[x[2], ]))
    }
    funCanberra <- function(x) {
        sum(abs(abund[x[1], ] - abund[x[2], col])/(abund[x[1], ] + abund[x[2], ]), na.rm=TRUE)
    }
    funBrayCurtis <- function(x) {
        sum(abs(abund[x[1], ] - abund[x[2], ])/sum(abund[c(x[1],x[2]), ]))
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    if (meth == "Euclidean")
        d <- unlist(apply(index, 1, funEuclidean))
    else if (meth == "Manhattan")
        d <- unlist(apply(index, 1, funManhattan))
    else if (meth == "Canberra")
        d <- unlist(apply(index, 1, funCanberra))
    else   
        d <- unlist(apply(index, 1, funBrayCurtis))
    attr(d, "Size") <- nrow(abund)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "method") <- meth
    class(d) <- "dist"
    return(d)
}

    distances <- distabund(comm, method)
    funSTATs <- function(dis, gro){
        QT <- sum(as.vector(dis))/Nplots
        Distances <- as.matrix(dis)
        ldistancesW <- lapply(levels(gro), function(x) Distances[gro==x, gro==x])
        ldistancesW <- lapply(ldistancesW, as.dist, diag = FALSE, upper = FALSE)
        lQW <- lapply(ldistancesW, function(x) sum(as.vector(x)/attributes(x)$Size))
        QW <- sum(as.vector(unlist(lQW)))
        QB <- QT - QW
        res <- c(QT, QB, QW)
        names(res) <- c("QT", "QB", "QW")
        return(res)
    }
    funperm <- function(i){
         egroups <- sample(groups)
         theo <- funSTATs(distances, egroups)[2]
         return(theo)
    }
    THEO <- sapply(1:nrep, funperm)
    OBS <- funSTATs(distances, groups)[2]
    simu <- as.randtest(sim = THEO, obs = OBS, alter = "greater")
    tabobs <- funSTATs(distances, groups)

}

else{

distabunde <- function (abund, col, meth) 
{

    d <- matrix(0, nrow(abund), nrow(abund))
    funEuclidean <- function(x) {
        (abund[x[1], col] - abund[x[2], col])^2
    }
    funManhattan <- function(x) {
        abs(abund[x[1], col] - abund[x[2], col])
    }
    funCanberra <- function(x) {
        if((abund[x[1], col] + abund[x[2], col]) < tol) return(0)
        else
        return(abs(abund[x[1], col] - abund[x[2], col])/(abund[x[1], col] + abund[x[2], col]))
    }
    funBrayCurtis <- function(x) {
        if(is.null(attributes(abund)$internal.p.rowSums))
        abs(abund[x[1], col] - abund[x[2], col])/sum(abund[c(x[1],x[2]), ])
        else
        abs(abund[x[1], col] - abund[x[2], col])/sum(attributes(abund)$rowSums[c(x[1],x[2])])
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    if (meth == "Euclidean")
        d <- unlist(apply(index, 1, funEuclidean))
    else if (meth == "Manhattan")
        d <- unlist(apply(index, 1, funManhattan))
    else if (meth == "Canberra")
        d <- unlist(apply(index, 1, funCanberra))
    else   
        d <- unlist(apply(index, 1, funBrayCurtis))
    attr(d, "Size") <- nrow(abund)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "method") <- meth
    class(d) <- "dist"
    return(d)
}

    listDIS <- lapply(1:ncol(comm), function(i) as.matrix(distabunde(comm, i, method)))
    if(global) {
        distances <- listDIS[[1]]
        for(i in 2:length(listDIS))    distances <- distances + listDIS[[i]]
    distances <- as.dist(distances, diag = FALSE, upper = FALSE)
    }
    listDIS <- lapply(listDIS, as.dist, diag = FALSE, upper = FALSE)
    funSTATs <- function(dis, gro){
        QT <- sum(as.vector(dis))/Nplots
        Distances <- as.matrix(dis)
        ldistancesW <- lapply(levels(gro), function(x) Distances[gro==x, gro==x])
        ldistancesW <- lapply(ldistancesW, as.dist, diag = FALSE, upper = FALSE)
        lQW <- lapply(ldistancesW, function(x) sum(as.vector(x)/attributes(x)$Size))
        QW <- sum(as.vector(unlist(lQW)))
        QB <- QT - QW
        res <- c(QT, QB, QW)
        names(res) <- c("QT", "QB", "QW")
        return(res)
    }
    funperm <- function(i){
         egroups <- sample(groups)
         if(global)   theoglobal <- funSTATs(distances, egroups)[2]
         theospecies <- as.vector(unlist(lapply(listDIS, function(dd) funSTATs(dd, egroups)[2])))
         if(global) theo <- c(theoglobal, theospecies)
         else theo <- theospecies
         return(theo)
    }
    THEO <- t(cbind.data.frame(sapply(1:nrep, funperm)))
    if(global) colnames(THEO) <- c("GLOBAL", colnames(comm))
    else colnames(THEO) <- colnames(comm)
    obsspecies <- as.vector(unlist(lapply(listDIS, function(dd) funSTATs(dd, groups)[2])))
    if(global){
    obsglobal <- funSTATs(distances, groups)[2]
    OBS <- c(obsglobal, obsspecies)
    }
    else OBS <- obsspecies
    simu <- as.krandtest(sim = THEO, obs = OBS, alter = "greater", p.adjust.method = padjust)
    obsspecies <- cbind.data.frame(lapply(listDIS, function(dd) funSTATs(dd, groups)))
    colnames(obsspecies) <- colnames(comm)
    if(global){
        obsglobal <- funSTATs(distances, groups)
        tabobs <- cbind.data.frame(Global = obsglobal, obsspecies)
    }
    else tabobs <- obsspecies
}
    RES <- list(observations = tabobs, test = simu, call = match.call(), method = method, padjust = padjust, tol = tol)
    class(RES) <- "dbMANOVAspecies"
    return(RES)
}
