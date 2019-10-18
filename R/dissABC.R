dissABC <-
function(comm, dis, option = 1:4, method = c("J", "S", "O", "K", "SS", "Si"))
{
    method <- method[1]
    if(!method%in%c("J", "S", "O", "K", "SS", "Si")) stop("method must be J, S, O, K, SS, or Si")
    if(inherits(dis, "dist")) dis <- as.matrix(dis)
    if(any(dis>1)){
        warning("Phylogenetic dissimilarities are not in the range 0-1. They have been normalized by the maximum")
        dis <- dis/max(dis)
    }
    if(any(!colnames(comm) %in%rownames(dis))) stop("At least one species in the matrix of abundances is missing in the matrix of dissimilarities") 
    if(any(!colnames(comm) %in%colnames(dis))) stop("At least one species in the matrix of abundances is missing in the matrix of dissimilarities") 
    dis <- dis[colnames(comm), colnames(comm)]
    dataset <- t(comm)
    similarities <- 1-as.matrix(dis)
    if(option[1]%in%c(1, 2)) rel_abu <- dataset
    else{
        total <- colSums(dataset)
        rel_abu <- sweep(dataset, 2, total, "/")
    }
    num.plot<-dim(dataset)[2]
    num.sp <- dim(dataset)[1]
    names<-list(colnames(dataset), colnames(dataset))
    sim.matrix<-matrix(0, nrow=num.plot, ncol=num.plot, dimnames=names)
    for (i in 2:num.plot) {
    for (j in 1:(i-1)) {
	    if(option[1]%in%c(1,3)) {	  
            garde <- (1:num.sp)[(rel_abu[, j]+rel_abu[, i])>0]
            sim <- similarities[garde, garde]
            x <- rel_abu[, j]
            x <- x[garde]
            y <- rel_abu[, i]
            y <- y[garde]
            mat_folk <- sim*x
            mat_folk2 <- sim*y
        }
	    else{
            mat_folk <- similarities*rel_abu[, j]
            mat_folk2 <- similarities*rel_abu[, i]
	    }
        Zik <- colSums(mat_folk)
        Zih <- colSums(mat_folk2)
        tabZ <- rbind.data.frame(Zik, Zih)
        A <- sum(sapply(tabZ, min))
        B <- sum(sapply(tabZ, max)-Zih)
        C <- sum(sapply(tabZ, max)-Zik)
        if(method == "J") index <- A/(A+B+C)
        else if(method == "S") index <- 2*A/(2*A+B+C)
        else if(method == "O") index <- A/(sqrt(A+B)*sqrt(A+C))
        else if(method == "K") index <- 0.5*(A/(A+B)+A/(A+C))
        else if(method == "SS") index <- A/(A+2*B+2*C)
        else if(method == "Si") index <- A/(A+min(c(B,C)))
        sim.matrix[i, j] <- index
    }
    }
    sim.matrix <- sim.matrix + t(sim.matrix)
    sim.matrix <- sim.matrix + diag(rep(1, num.plot))
    return(sim.matrix)
}
