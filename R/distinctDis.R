distinctDis <-
function(dis, method=c("Rb","AV","FV","NN","full"), standardized = FALSE){

    if (any(is.na(match(method, c("Rb","AV","FV","NN","full"))))) 
        stop("unconvenient method")
    if(any(method=="full")) method <- c("Rb","AV","FV","NN")
    nbMeth <- length(method)
    dis <- as.matrix(dis)
    nsp <- nrow(dis)
    resWeights <- as.data.frame(matrix(0, nsp, nbMeth))
    rownames(resWeights) <- rownames(dis)
    for (k in 1:nbMeth) {
        meth <- method[k]
        if (meth == "Rb") {
            ori <- abs(eigen(dis)$vector[, 1])
            if(standardized)
                ori <- ori/sum(ori)
            resWeights[, k] <- ori
            names(resWeights)[k] <- "Rb"
        }
        if (meth == "AV") {
            fun <- function(i){
                return(mean(dis[i, -i]))
            }
            ori <- sapply(1:nsp, fun)
            if(standardized)
                ori <- ori/sum(ori)
            resWeights[, k] <-  ori
            names(resWeights)[k] <- "AV"
        }
        if (meth == "FV") {
            fun <- function(i){
                return(mean(dis[i, ]))
            }
            ori <- sapply(1:nsp, fun)
            if(standardized)
                ori <- ori/sum(ori)
            resWeights[, k] <-  ori
            names(resWeights)[k] <- "FV"
        }
        if (meth == "NN") {
            fun <- function(i){
                return(min(dis[i, -i]))
            }
            ori <- sapply(1:nsp, fun)
            if(standardized)
                ori <- ori/sum(ori)
            resWeights[, k] <-  ori
            names(resWeights)[k] <- "NN"
        }

    }
    return(resWeights)
    
}
