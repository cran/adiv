eveparam <-
function(comm, method = c("hill", "tsallis","renyi"), q = 2, tol = 1e-8){

  m <- comm
  method <- method[1]
  if(!method%in%c("tsallis", "hill", "renyi")) stop("Incorrect method")
  nsp <- ncol(m)
  ncom <- nrow(m)
  if(any(m<0)) stop("m should contain nonnegative values")
  if(any(rowSums(m)==0)) stop("empty communities should be discarded")
  m <- sweep(m, 1, rowSums(m), "/")
  m <- t(m)

  tsallis <- function(x, q){
        funtsallis <- function(y, q) {
            y <- y[y>0]
            if(length(y)==1) return(0)
            if(abs(q-1) < tol){
                resi <- -sum(y*log(y))/log(length(y))
            }
            else{
            	resi <- (1-sum(y^q))/(1-(length(y))^(1-q))
            }
            return(resi)
        }
        if(ncol(x)==1) res <- funtsallis(x[, 1], q)
        else res <- apply(x, 2, funtsallis, q)
        return(res)
    }
    hill <- function(x, q){
        funhill <- function(y, q) {
            y <- y[y>0]
            if(abs(q-1) < tol){
                resi <- exp(-sum(y*log(y)))/length(y)
            }
            else{
            	resi <- (sum(y^q))^(1/(1-q))/length(y)
            }
            return(resi)
        }
        if(ncol(x)==1) res <- funhill(x[, 1], q)
        else res <- apply(x, 2, funhill, q)
        return(res)

    }
    renyi <- function(x, q){
        funrenyi <- function(y, q) {
            y <- y[y>0]
            if(abs(q-1) < tol){
                resi <- -sum(y*log(y))/log(length(y))
            }
            else{
            	resi <- log((sum(y^q))^(1/(1-q)))/log(length(y))
            }
            return(resi)
        }
        if(ncol(x)==1) res <- funrenyi(x[, 1], q)
        else res <- apply(x, 2, funrenyi, q)
        return(res)
    }
    if( length(q)==1 ){
        if(method == "tsallis"){
            vres <- tsallis(m, q)
            class(vres) <- "eveparam"
            return(vres)
        }
        if(method == "hill"){
            vres <- hill(m, q) 
            class(vres) <- "eveparam"
            return(vres)
        }
        if(method == "renyi"){
            vres <- renyi(m, q) 
            class(vres) <- "eveparam"
            return(vres)
        }
    }
    if ( length(q) > 1){
        if(method == "tsallis"){
           calcul1 <- sapply(q, function(x) tsallis(m, x))
           if(ncom>1)
               tab1 <- cbind.data.frame(calcul1)
           else
               tab1 <- as.data.frame(matrix(calcul1, 1, byrow=TRUE))
           rownames(tab1) <- colnames(m)
           listtotale <- list()
           listtotale$q <- q
           listtotale$eve <- tab1
           class(listtotale) <- "eveparam"
           return(listtotale)
        }
        if(method == "hill"){
           calcul1 <- sapply(q, function(x) hill(m, x))
           if(ncom>1)
               tab1 <- cbind.data.frame(calcul1)
           else
               tab1 <- as.data.frame(matrix(calcul1, 1, byrow=TRUE))
           rownames(tab1) <- colnames(m)
           listtotale <- list()
           listtotale$q <- q
           listtotale$eve <- tab1
           class(listtotale) <- "eveparam"
           return(listtotale)
        }
        if(method == "renyi"){
           calcul1 <- sapply(q, function(x) renyi(m, x))
           if(ncom>1)
               tab1 <- cbind.data.frame(calcul1)
           else
               tab1 <- as.data.frame(matrix(calcul1, 1, byrow=TRUE))
           rownames(tab1) <- colnames(m)
           listtotale <- list()
           listtotale$q <- q
           listtotale$eve <- tab1
           class(listtotale) <- "eveparam"
           return(listtotale)
        }
    }
  

}
