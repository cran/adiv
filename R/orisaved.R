orisaved <-
function(phyl, rate = 0.1, method = 1:2)
{
  arg.phyl <- .checkphyloarg(phyl)
  phyl <- arg.phyl$phyl
  phyl.phylo <- arg.phyl$phyl.phylo
  rm(arg.phyl)
  
  method <- method[1]
  if (any(is.na(match(method, 1:2)))) stop("unconvenient method")
  if (length(rate) != 1) stop("unconvenient rate")
  if (!is.numeric(rate)) stop("rate must be a real value")
  if (!(rate>=0 & rate<=1)) stop("rate must be between 0 and 1")
  if (rate == 0) return(0)

  if(is.binary.phylo(phyl.phylo))
      phy.h <- as.hclust(phyl.phylo) ## also test for ultrametricity 
  else{
      if(!is.ultrametric(phyl.phylo)) stop("the tree is not ultrametric")
      phy.h <- hclust(as.dist(cophenetic.phylo(phyl.phylo)/2), "average")
  }
  phyl.D <- cophenetic.phylo(phyl.phylo)/2
  
  nbesp <- nTips(phyl)
  
  Rate <- round(seq(0, nbesp, by = nbesp * rate))
  Rate <- Rate[-1]
  num.Orig <- as.vector(solve(phyl.D, rep(1, nbesp)))
  denum.Orig <- as.vector(t(rep(1, nbesp))%*%num.Orig)
  Orig <- as.vector(num.Orig/denum.Orig)
  
  OrigCalc <- function(i) {
    if (method == 1) {
      return(sum(unlist(lapply(split(Orig, cutree(phy.h, i)), max))))
    }
    if (method == 2) {
      return(sum(unlist(lapply(split(Orig, cutree(phy.h, i)), min))))
    }
  }
  res <- c(0, sapply(Rate, OrigCalc))
  return(res)
}
