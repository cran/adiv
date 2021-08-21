distinctTree <- function (phyl, method = c("ED", "ES", "Delta*"), palpha = 0, standardized = FALSE) 
{
    arg.phyl <- .checkphyloarg(phyl)
    phyl.phylo <- arg.phyl$phyl.phylo
    phyl <- as(collapse.singles(phyl.phylo, root.edge=TRUE), "phylo4")
    rm(arg.phyl)
    if (any(!method %in% c("ED", "ES", "Delta*"))) 
        stop("unconvenient method")
    nbMeth <- length(method)
    nbesp <- nTips(phyl)
    nbnodes <- nNodes(phyl)
    if("Delta"%in%method) nbMeth <- nbMeth + length(palpha)-1
    resWeights <- as.data.frame(matrix(0, nbesp, nbMeth))
    rownames(resWeights) <- tipLabels(phyl)
    for (k in 1:length(method)) {
        meth <- method[k]
        if (meth == "ED") {
            allPath <- lapply(getNode(phyl, type = "tip"), 
                function(tip) .tipToRoot(phyl, tip, getNode(phyl, 
                  nbesp + 1), include.root = FALSE))
            nTip.node <- lapply(listTips(phyl), length)
            length.node <- edgeLength(phyl, getNode(phyl, type = "internal"))
            weight.node <- length.node/unlist(nTip.node)
            reslist <- lapply(allPath, function(x) sum(unlist(weight.node)[x - 
                nbesp]))
            res <- unlist(reslist) + edgeLength(phyl, getNode(phyl, 
                type = "tip"))
            if (standardized) 
                res <- res/sum(res)
            resWeights[, k] <- res
            names(resWeights)[k] <- "ED"
        }
        else if(meth == "ES"){
            allPath <- lapply(getNode(phyl, type = "tip"), 
                function(tip) .tipToRoot(phyl, tip, getNode(phyl, 
                  nbesp + 1), include.root = FALSE))
            ndescendents <- sapply(listDD(phyl), length)
            length.node <- edgeLength(phyl, getNode(phyl, type = "internal"))
            reslist <- lapply(allPath, function(x) sum(length.node[x - 
                nbesp]/cumprod(ndescendents[x - nbesp]), na.rm = TRUE))
            res <- unlist(reslist) + edgeLength(phyl, getNode(phyl, 
                type = "tip"))
            if (standardized) 
                res <- res/sum(res)
            resWeights[, k] <- res
            names(resWeights)[k] <- "ES"
        }
        else{
            allPath <- lapply(getNode(phyl, type = "tip"), 
                function(tip) .tipToRoot(phyl, tip, getNode(phyl, 
                  nbesp + 1), include.root = FALSE))
            nTip.node <- unlist(lapply(listTips(phyl), length))
            length.node <- edgeLength(phyl, getNode(phyl, type = "internal"))
            if(length(palpha)<2){
                if(abs(palpha-1)>1e-10)
                reslist <- lapply(allPath, function(x) sum(length.node[x - 
                nbesp]*(nbesp^(palpha-1) - nTip.node[x - nbesp]^(palpha-1))/(nbesp^(palpha-1) - 1), na.rm = TRUE))
                else
                reslist <- lapply(allPath, function(x) sum(length.node[x - 
                nbesp]*(1-log(nTip.node[x - nbesp], nbesp)), na.rm = TRUE))
                res <- unlist(reslist) + edgeLength(phyl, getNode(phyl, 
                type = "tip"))
                if (standardized) 
                res <- res/sum(res)
                resWeights[, k] <- res
                names(resWeights)[k] <- "Delta"
            }
            else{
                for(i in 1:length(palpha)){
                if(abs(palpha[i]-1)>1e-10)
                reslist <- lapply(allPath, function(x) sum(length.node[x - 
                nbesp]*(nbesp^(palpha[i]-1) - nTip.node[x - nbesp]^(palpha[i]-1))/(nbesp^(palpha[i]-1) - 1), na.rm = TRUE))
                else
                reslist <- lapply(allPath, function(x) sum(length.node[x - 
                nbesp]*(1-log(nTip.node[x - nbesp], nbesp)), na.rm = TRUE))
                res <- unlist(reslist) + edgeLength(phyl, getNode(phyl, 
                type = "tip"))
                if (standardized) 
                res <- res/sum(res)
                resWeights[, k+i-1] <- res
                names(resWeights)[k+i-1] <- paste("Delta", i, sep="_")
                }
            }
        }
    }
    return(resWeights)
}
