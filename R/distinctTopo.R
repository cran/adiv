distinctTopo <-
function (phyl, method = 1:3, standardized = FALSE) 
{
  arg.phyl <- .checkphyloarg(phyl)
  phyl <- arg.phyl$phyl
  phyl.phylo <- arg.phyl$phyl.phylo
  rm(arg.phyl)

  ## phyl is a phylo4 object, phyl.phylo is a phylo object

    if (any(is.na(match(method, 1:3)))) 
        stop("unconvenient method")
    nbMeth <- length(method)
    nbesp <- nTips(phyl)
    nbnodes <- nNodes(phyl)
    resWeights <- as.data.frame(matrix(0, nbesp, nbMeth))
    rownames(resWeights) <- tipLabels(phyl)
    for (k in 1:nbMeth) {
        meth <- method[k]
        if (meth == 1) {
          interm <- distRoot(phyl, method = "nNodes")
          res <- 1/interm/min(1/interm)
          if(standardized)
                res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "VW"
        }
        if (meth == 2) {
          interm <- distRoot(phyl, method = "sumDD")
          res <- 1/interm/min(1/interm)
          if(standardized)
                res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "M"
        }
        if (meth == 3) {
          interm <- distRoot(phyl, method = "Abouheif")
          res <- 1/interm/min(1/interm)
          if(standardized)
                res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "A"          
        }
      }
  return(resWeights)
}
