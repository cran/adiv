distinctTree <-
function(phyl, method=c("ED","ES"), standardized = FALSE) 
{
  arg.phyl <- .checkphyloarg(phyl)
  phyl <- arg.phyl$phyl
  phyl.phylo <- arg.phyl$phyl.phylo
  rm(arg.phyl)

  ## phyl is a phylo4 object, phyl.phylo is a phylo object

    if (any(is.na(match(method, c("ED","ES"))))) 
        stop("unconvenient method")
    nbMeth <- length(method)
    nbesp <- nTips(phyl)
    nbnodes <- nNodes(phyl)
    resWeights <- as.data.frame(matrix(0, nbesp, nbMeth))
    rownames(resWeights) <- tipLabels(phyl)
    for (k in 1:nbMeth) {
        meth <- method[k]
        if (meth == "ED") {
    allPath <- lapply(getNode(phyl, type = "tip"), function(tip) .tipToRoot(phyl, tip, getNode(phyl, nbesp + 1), include.root = FALSE))
    nTip.node <- lapply(listTips(phyl), length)
    length.node <- edgeLength(phyl, getNode(phyl, type = "internal"))
    weight.node <- length.node / unlist(nTip.node)
    reslist <- lapply(allPath,function(x) sum(unlist(weight.node)[x - nbesp]))
    res <- unlist(reslist) + edgeLength(phyl, getNode(phyl, type = "tip"))
    if(standardized)
         res <- res/sum(res)
    resWeights[, k] <- res
    names(resWeights)[k] <- "ED"
        }
        else{
          allPath <- lapply(getNode(phyl, type = "tip"), function(tip) .tipToRoot(phyl, tip, getNode(phyl, nbesp + 1), include.root = FALSE))
          ndescendents <- sapply(listDD(phyl),length)
          length.node <- edgeLength(phyl, getNode(phyl, type = "internal"))
          reslist <- lapply(allPath, function(x) 
                            sum(length.node[x - nbesp] / cumprod(ndescendents[x - nbesp]), na.rm = TRUE))
          res <- unlist(reslist) + edgeLength(phyl, getNode(phyl, type = "tip"))
          if(standardized)
              res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "ES"
        }
     }
     return(resWeights)
}
