uniqueness <-
function(comm, dis, tol = 1e-8, abundance = TRUE){

	if(!is.null(colnames(comm)) & !is.null(attributes(dis)$Labels)){
		if(any(!colnames(comm)%in%attributes(dis)$Labels)) stop("One or several species in comm are not in dis; check species names in comm and in dis")
		else dis <- as.dist(as.matrix(dis)[colnames(comm), colnames(comm)])
	}
	else if(ncol(comm)!=attributes(dis)$Size) stop("the number of species in comm must be equal to that in dis") 		

	D <- as.matrix(dis)
     if(any(D>1)){
         D <- D/max(D)
         warnings("All values in dis have been divided by their maximum (highest observed value in D)")
     }

	if(!abundance){
		comm[comm>0] <- 1
	}
        commt <- as.data.frame(t(comm))

        funK <- function(v){
		p <- v/sum(v)
		K <- apply(D, 1, function(x) sum(x*p))
		K[p<tol] <- NA
                return(K)
        }
        V <- cbind.data.frame(sapply(commt, funK))
	rownames(V) <- colnames(comm)
	colnames(V) <- rownames(comm)
        funKbar <- function(v){
		p <- v/sum(v)
		Kbar <- sapply(1:nrow(D), function(i) sum(D[i,]*v/sum(v[-i])))
 		Kbar[p<tol] <- NA
                return(Kbar)
        }
	Kbar <- cbind.data.frame(sapply(commt, funKbar))
	rownames(Kbar) <- colnames(comm)
	colnames(Kbar) <- rownames(comm)
        funQ <- function(v){
		p <- v/sum(v)
                Q <- t(p)%*%D%*%p
		return(Q)
	}
        Q <- unlist(sapply(commt, funQ))

        funSim <- function(v){
		p <- v/sum(v)
                S <- 1-sum(p^2)
                return(S)
	}
        Sim <- unlist(sapply(commt, funSim))

        funN <- function(v){
		p <- v/sum(v)
                N <- length(p[p>tol])
                return(N)
	}
        N <- unlist(sapply(commt, funN))
	   U <- Q/Sim
        R <- 1- U
        Ustar <- (1-Sim)/(1-Q)
        Rstar <- 1-Ustar


	red <- cbind.data.frame(N=N, D=Sim, Q=Q, U=U, R=R, Ustar=Ustar, Rstar=Rstar)
	rownames(red) <- rownames(comm)

        res <- list()
        res$Kbar <- Kbar
        res$V <- V
	res$red <- red
        return(res)

}
