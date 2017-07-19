plot.aptree <-
function (x, col.line = 'blue', ...)
{

  if(!inherits(x, "aptree"))
    stop("x must be an object of class aptree")

  thecall <- as.list(attributes(x)$call)
  tre <- eval(thecall$phyl, sys.frame(0))
  tre <- .checkphyloarg(tre)
  phyl.phylo <- tre$phyl.phylo
  com <- eval(thecall$comm, sys.frame(0))
  if(length(phyl.phylo$tip.label)>ncol(com))
      phyl.phylo <- drop.tip(phyl.phylo, phyl.phylo$tip.label[!phyl.phylo$tip.label%in%colnames(com)])    
    
  plot(phyl.phylo, type = "phylogram", ...)
  theplot <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  ord <- as.numeric(substr(x$names, 2, nchar(x$names)))
  xx <- sort(theplot$xx, decreasing = TRUE)
  xx <- xx[-(1:(nTips(phyl.phylo)+1))]
  abline(v = xx, col = col.line, lty = 2)
  
}
