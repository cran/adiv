distinctTopo <-
function (phyl, method = c("VW","M","A","full"), standardized = FALSE) 
{
  arg.phyl <- .checkphyloarg(phyl)
  phyl <- arg.phyl$phyl
  phyl.phylo <- arg.phyl$phyl.phylo
  rm(arg.phyl)

  ## phyl is a phylo4 object, phyl.phylo is a phylo object

    if (any(is.na(match(method, c("VW","M","A","full"))))) 
        stop("unconvenient method")
    if(any(method=="full")) method <- c("VW","M","A")
    nbMeth <- length(method)
    nbesp <- nTips(phyl)
    nbnodes <- nNodes(phyl)
    resWeights <- as.data.frame(matrix(0, nbesp, nbMeth))
    rownames(resWeights) <- tipLabels(phyl)
    for (k in 1:nbMeth) {
        meth <- method[k]
        if (meth == "VW") {
          interm <- distRoot(phyl, method = "nNodes")
          res <- 1/interm/min(1/interm)
          if(standardized)
                res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "VW"
        }
        if (meth == "M") {
          interm <- distRoot(phyl, method = "sumDD")
          res <- 1/interm/min(1/interm)
          if(standardized)
                res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "M"
        }
        if (meth == "A") {
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
