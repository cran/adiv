distinctDis <- function (dis, method = c("Rb", "AV", "FV", 
    "NN", "Dstar", "full"), palpha = 0, standardized = FALSE) 
{
    if (any(is.na(match(method, c("Rb", "AV", "FV", 
        "NN", "Dstar", "full"))))) 
        stop("unconvenient method")
    if (any(method == "full")) 
        method <- c("Rb", "AV", "FV", "NN", "Dstar")
    nbMeth <- length(method)
    dis <- as.matrix(dis)
    nsp <- nrow(dis)
    if("D"%in%method) nbMeth <- nbMeth + length(palpha)-1
    if("Dstar"%in%method) nbMeth <- nbMeth + length(palpha)-1
    resWeights <- as.data.frame(matrix(0, nsp, nbMeth))
    rownames(resWeights) <- rownames(dis)
    for (k in 1:length(method)) {
        meth <- method[k]
        if (meth == "Rb") {
            ori <- abs(eigen(dis)$vector[, 1])
            if (standardized) 
                ori <- ori/sum(ori)
            resWeights[, k] <- ori
            names(resWeights)[k] <- "Rb"
        }
        if (meth == "AV") {
            fun <- function(i) {
                return(mean(dis[i, -i]))
            }
            ori <- sapply(1:nsp, fun)
            if (standardized) 
                ori <- ori/sum(ori)
            resWeights[, k] <- ori
            names(resWeights)[k] <- "AV"
        }
        if (meth == "FV") {
            fun <- function(i) {
                return(mean(dis[i, ]))
            }
            ori <- sapply(1:nsp, fun)
            if (standardized) 
                ori <- ori/sum(ori)
            resWeights[, k] <- ori
            names(resWeights)[k] <- "FV"
        }
        if (meth == "NN") {
            fun <- function(i) {
                return(min(dis[i, -i]))
            }
            ori <- sapply(1:nsp, fun)
            if (standardized) 
                ori <- ori/sum(ori)
            resWeights[, k] <- ori
            names(resWeights)[k] <- "NN"
        }
        if (meth == "Dstar") {
            dissort <- apply(dis, 1, sort)
            dissort <- t(apply(dissort, 2, diff))
            funDstar <- function(j, pa) {
                if(abs(pa-1)<1e-10)
                return( sum( dissort[j, ]*(1-log(1:(nsp-1),nsp)) ) )
                else
                return(sum(dissort[j, ]*(nsp^(pa-1)-(1:(nsp-1))^(pa-1))/(nsp^(pa-1)-1)))
            }
            if(length(palpha)<2){
                ori <- sapply(1:nsp, function(x) funDstar(x, palpha))
                if (standardized) 
                ori <- ori/sum(ori)
                resWeights[, k] <- ori
                names(resWeights)[k] <- "Dstar"
            }
            else {
               for(i in 1:length(palpha)){
                ori <- sapply(1:nsp, function(x) funDstar(x, palpha[i]))
                if (standardized) 
                ori <- ori/sum(ori)
                resWeights[, k+i-1] <- ori
                names(resWeights)[k+i-1] <- paste("Dstar", i, sep="_")
               }
            }
        }

    }
    return(resWeights)
}
