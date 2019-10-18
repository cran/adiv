distinctUltra <-
function (phyl, method = c("Qb","2Hb")) 
{
  arg.phyl <- .checkphyloarg(phyl)
  phyl.phylo <- arg.phyl$phyl.phylo
  phyl <- arg.phyl$phyl
  ## phyl is a phylo4 object, phyl.phylo is a phylo object

   if (any(is.na(match(method, c("Qb","2Hb"))))) 
        stop("unconvenient method")
    nbMeth <- length(method)
    nbesp <- nTips(phyl)
    nbnodes <- nNodes(phyl)
    resWeights <- as.data.frame(matrix(0, nbesp, nbMeth))
    rownames(resWeights) <- tipLabels(phyl)
    for (k in 1:nbMeth) {
        meth <- method[k]
        if (meth == "Qb") {
          if(!is.ultrametric(phyl.phylo))
            stop("phyl must be an ultrametric tree")
          D <- cophenetic.phylo(phyl.phylo)/2
          num.Orig <- as.vector(solve(D, rep(1, nbesp)))
          denum.Orig <- as.vector(t(rep(1, nbesp))%*%num.Orig)
          res <- num.Orig/denum.Orig
          resWeights[, k] <- res
          names(resWeights)[k] <- "QEbased"
        }
        if (meth == "2Hb") {
          if(!is.ultrametric(phyl.phylo))
            stop("phyl must be an ultrametric tree")
          C <- vcv.phylo(phyl.phylo)
          C2 <- C^2
          d <- diag(C)
          num.Orig <- (solve(C2)%*%d)
          denum.Orig <- sum(num.Orig)
          res <- num.Orig/denum.Orig
          resWeights[, k] <- res
          names(resWeights)[k] <- "twoHbased"
        }
     }
     return(resWeights)
}
