betastatsor <-
function(comm){ 
    comm <- ifelse(comm > 0, 1, 0) 
    N <- nrow(comm) 
    S <- ncol(comm[, colSums(comm)>0]) 
    sumNi <- sum(comm) 
    R <- rowSums(comm) 
    spmin <- apply(comm, 2, function(x) min((x * R)[x > 0])) 
    n0 <- spmin 
    for (i in 1:ncol(comm)) n0[i] <- sum(comm[, i] == 0 & R >= 
        spmin[i]) 
    NiT <- sum(n0) 
    beta <- (S*N-sumNi)/(sumNi*(N-1)) 
    sim <- 1-beta 
    betaT <- NiT/(sumNi*(N-1)) 
    betaN <- beta - betaT 
    res <- c(beta, betaT, betaN, sim) 
    names(res) <- c("beta", "betaT", "betaN", "sim") 
    return(res) 
}
