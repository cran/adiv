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
    CtipToRootIncl <- function(TRE4, TIPS){
        LL <- lapply(TIPS, function(ti) ancestors(TRE4, ti))
        return(LL)
    }
    dissRoot <- function(tre4, method = c("nNodes", "sumnChildren", "prodnChildren")){
        allPath <- CtipToRootIncl(tre4, getNode(tre4, type = "tip"))
        if(method == "nNodes") { 
        res <- as.vector(unlist(lapply(allPath, length)))
        names(res) <- names(allPath)
        return(res)
        }
        if(method == "sumnChildren") { 
        FUNchildren <- function(No){
            return(sapply(No, function(Noi) length(descendants(tre4, Noi, "children"))))
        } 
        LL <- lapply(allPath, FUNchildren)
        res <- as.vector(unlist(lapply(LL, sum)))
        names(res) <- names(allPath)
        return(res)
        }
        if(method == "prodnChildren") { 
        FUNchildren <- function(No){
            return(sapply(No, function(Noi) length(descendants(tre4, Noi, "children"))))
        } 
        LL <- lapply(allPath, FUNchildren)
        res <- as.vector(unlist(lapply(LL, prod)))
        names(res) <- names(allPath)
        return(res)
        }

    }


    for (k in 1:nbMeth) {
        meth <- method[k]
        if (meth == "VW") {
          interm <- dissRoot(phyl, method = "nNodes")
          res <- 1/interm/min(1/interm)
          if(standardized)
                res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "VW"
        }
        if (meth == "M") {
          interm <- dissRoot(phyl, method = "sumnChildren")
          res <- 1/interm/min(1/interm)
          if(standardized)
                res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "M"
        }
        if (meth == "A") {
          interm <- dissRoot(phyl, method = "prodnChildren")
          res <- 1/interm/min(1/interm)
          if(standardized)
                res <- res/sum(res)
          resWeights[, k] <- res
          names(resWeights)[k] <- "A"          
        }
      }
  return(resWeights)
}
