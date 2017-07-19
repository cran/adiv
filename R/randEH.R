randEH <-
function(phyl, nbofsp, nrep = 10)
{
  if (length(nbofsp)!= 1) stop("unconvenient nbofsp")
  arg.phyl <- .checkphyloarg(phyl)
  phyl <- arg.phyl$phyl
  phyl.phylo <- arg.phyl$phyl.phylo
  rm(arg.phyl)
  
  nbesp <- nTips(phyl)
  if (!((0 <= nbofsp) & (nbofsp <= nbesp)))
    stop("unconvenient nbofsp")
  nbofsp <- round(nbofsp)
  if (nbofsp == 0) return(rep(0, nrep))
  if (nbofsp == nbesp) {
    return(rep(sumEdgeLength(phyl), nrep))
  }
  phyl.D <- NULL
  if(is.ultrametric(phyl.phylo)) {
    phyl.D <- cophenetic.phylo(phyl.phylo)/2
    
    EvoHist.util <- function(phyl, select, phylD) {
      if(length(select)==1)
        sum1 <- max(phyl.D)
      else if(length(select)==2) {
        sum1 <- phyl.D[select[1], select[2]] + max(phyl.D)
      } else {
        fun.EvoHist <- function(i) {
          min(phyl.D[select[i], select[1:(i - 1)]])
        }
        sum1 <-  phyl.D[select[1], select[2]] + max(phyl.D) + sum(sapply(3:length(select), fun.EvoHist))
      }
      return(sum1)
    }
  } else {
    
    EvoHist.util <- function(phyl, select, phylD){
      sum1 <- sumEdgeLength(phyl, unique(unlist(ancestors(phyl, select, type="ALL"))))
      return(sum1)
    }
  }
  
  simuA1 <- function(i, phy, phylD) {
    comp <- sample(1:nbesp, nbofsp)
    resc <- EvoHist.util(phyl, select = comp, phylD)
    return(resc)
  }
  res <- sapply(1:nrep, simuA1, phyl, phyl.D)
  return(res)
}
