plot.rtestaptree <-
function (x, col.line = c('blue', 'red'), alpha = 0.05, ...)
{

  if(!inherits(x, "rtestaptree"))
    stop("x must be an object of class rtestaptree")

  thecall <- as.list(x$call)
  tre <- eval(thecall$phyl, sys.frame(0))
  tre <- .checkphyloarg(tre)
  phyl.phylo <- tre$phyl.phylo
  com <- eval(thecall$comm, sys.frame(0))
  if(length(phyl.phylo$tip.label)>ncol(com))
      phyl.phylo <- drop.tip(phyl.phylo, phyl.phylo$tip.label[!phyl.phylo$tip.label%in%colnames(com)])    
    
  plot(phyl.phylo, type = "phylogram", ...)
  theplot <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  ord <- as.numeric(substr(x$names, 2, nchar(x$names)))
  idx <- x$pvalue<=alpha
  xx <- sort(theplot$xx, decreasing = TRUE)
  xx <- xx[-(1:(nTips(phyl.phylo)+1))]
  abline(v = xx, col = ifelse(idx, col.line[2], col.line[1]), lty = 2)
  
}
