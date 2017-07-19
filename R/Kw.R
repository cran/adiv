Kw <-
function(phyl, trait, nrep = 999, alter = c("greater", "less", "two-sided")){

    alter <- alter[1]
    arg.phyl <- .checkphyloarg(phyl)
    phy <- arg.phyl$phyl.phylo
    rm(arg.phyl)

    mat <- vcv.phylo(phy)
    vars <- diag(mat)
    ntax <- length(phy$tip.label)
    Cm <- solve(mat)
    diagmat2 <- diag((diag(as.vector(t(rep(1, ntax))%*%Cm))%*%mat)/sum(Cm))

    KcalcIw <- function(x) 
    {
	
	ahat <- as.vector((t(rep(1, ntax))%*%Cm%*%x)/sum(Cm))
	xa <-x-ahat
	MSE0.MSEobs <- ((t(rep(1, ntax))%*%Cm%*%xa^2)/sum(Cm))/(t(xa)%*%Cm%*%xa)
  	MSE0.MSE <- 1/(ntax - 1) * (sum(diagmat2) - 1/sum(Cm))
  	K <- MSE0.MSEobs/MSE0.MSE
  	return(K)
    }

    Fun0 <- function(i){
	return(as.vector(KcalcIw(sample(trait))))
    }
    obs <- as.vector(KcalcIw(trait))
    executFun0 <- sapply(1:nrep, Fun0)
    res <- as.randtest(obs=obs, sim=executFun0, alter=alter)
    return(res)

}
