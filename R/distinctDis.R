distinctDis <-
function(dis, method=1:3, standardized = FALSE){

    if (any(is.na(match(method, 1:3)))) 
        stop("unconvenient method")
    nbMeth <- length(method)
    dis <- as.matrix(dis)
    nsp <- nrow(dis)
    resWeights <- as.data.frame(matrix(0, nsp, nbMeth))
    rownames(resWeights) <- attributes(dis)$Labels
    for (k in 1:nbMeth) {
        meth <- method[k]
        if (meth == 1) {
            ori <- abs(eigen(dis)$vector[, 1])
            if(standardized)
                ori <- ori/sum(ori)
            resWeights[, k] <- ori
            names(resWeights)[k] <- "Rb"
        }
        if (meth == 2) {
            fun <- function(i){
                return(mean(dis[i, -i]))
            }
            ori <- sapply(1:nsp, fun)
            if(standardized)
                ori <- ori/sum(ori)
            resWeights[, k] <-  ori
            names(resWeights)[k] <- "AV"
        }
        if (meth == 3) {
            fun <- function(i){
                return(mean(dis[i, ]))
            }
            ori <- sapply(1:nsp, fun)
            if(standardized)
                ori <- ori/sum(ori)
            resWeights[, k] <-  ori
            names(resWeights)[k] <- "FV"
        }
    }
    return(resWeights)
    
}
