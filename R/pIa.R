pIa <-
function(phyl, comm, exponent = 2, tol = 1e-8){

  if(is.vector(comm))
    return(sum(aptree(phyl, comm, exponent, tol)))
  else{
    res <- colSums(aptree(phyl, comm, exponent, tol))
    tab <- cbind.data.frame(res)
    names(tab) <- "diversity"
    rownames(tab) <- rownames(comm)
    return(tab)
  }
  
}
