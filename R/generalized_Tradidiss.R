generalized_Tradidiss <-
function(comm, dis, method = c("GC", "MS", "PE"), abundance = c("relative", "absolute", "none"), weights = c("uneven", "even"), tol = 1e-8)
{
    if(inherits(dis, "dist")) dis <- as.matrix(dis)
    if(any(dis>1)){
        warning("Phylogenetic dissimilarities are not in the range 0-1. They have been normalized by the maximum")
        dis <- dis/max(dis)
    }
    if(any(!colnames(comm) %in%rownames(dis))) stop("At least one species in the matrix of abundances is missing in the matrix of dissimilarities")
    if(any(!colnames(comm) %in%colnames(dis))) stop("At least one species in the matrix of abundances is missing in the matrix of dissimilarities")
    method <- method[1]
    if(!method%in%c("GC", "MS", "PE")) stop("Incorrect definition for method")
    abundance <- abundance[1]
    if(!abundance%in%c("relative", "absolute", "none"))
stop("Incorrect definition for abundance")
    weights <- weights[1]
    if(!weights%in%c("uneven", "even")) stop("Incorrect definition for weights")
dis <- dis[colnames(comm), colnames(comm)]
    dataset <- t(comm)
    similarities <- 1-as.matrix(dis)
total <- colSums(dataset)
    if(abundance == "relative")
        abu <- sweep(dataset, 2, total, "/")
    else if(abundance == "none"){
abu <- dataset
        abu[abu>tol] <- 1
        abu[abu<=tol] <- 0
}
    else abu <- dataset
    num.plot<-dim(dataset)[2]
    num.sp <- dim(dataset)[1]
    names<-list(colnames(dataset), colnames(dataset))
    dis.matrix<-matrix(0, nrow=num.plot, ncol=num.plot, dimnames=names)
    for (i in 2:num.plot) {
    for (j in 1:(i-1)) {
        mat_folk <- similarities*abu[, j]
        mat_folk2 <- similarities*abu[, i]
        Zik <- colSums(mat_folk)
        Zih <- colSums(mat_folk2)
        tabZ <- rbind.data.frame(Zik, Zih)
        garde <- apply(tabZ, 2, sum)>tol & apply(abu[,c(i,j)], 1, sum)>tol
Zik <- Zik[garde]
        Zih <- Zih[garde]
tabZ <- tabZ[, garde]
        if(weights == "even") wk <- rep(1/length(Zih), length(Zih))
        else {
        wk <- (abu[, j]+abu[, i])/sum(abu[, c(i,j)])
wk <- wk[garde]
        }
        if(method == "GC") index <- sum((wk*abs(Zik-Zih)/apply(tabZ, 2, sum)))
        else if(method == "MS") index <- sum((wk*abs(Zik-Zih)/apply(tabZ, 2, max)))
        else {
            zik <- Zik / (Zik + Zih)
            zih <- Zih / (Zik + Zih)
            tabz <- rbind.data.frame(zik, zih)
            enlev <- as.vector(apply(tabz, 2,
                function(x) any(x<=tol)))
            zik <- zik[!enlev]
            zih <- zih[!enlev]
            wkN <- wk[!enlev]
            wkY <- wk[enlev]
index <- sum(wkY) + sum(wkN*(1
                +(zik * log(zik, base = 2)
                + zih * log(zih, base = 2))))
        }
        dis.matrix[i, j] <- index
}
    }
    dis.matrix <- dis.matrix + t(dis.matrix)
    dis.matrix <- dis.matrix + diag(rep(0, num.plot))
return(as.dist(dis.matrix))
}
