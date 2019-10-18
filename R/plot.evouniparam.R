plot.evouniparam <-
function(x, legend = TRUE, legendposi = "topright", axisLABEL = "Tree-based uniqueness", type="b", col = if(is.numeric(x)) NULL else sample(colors(distinct = TRUE), nrow(x$uni)), lty = if(is.numeric(x)) NULL else rep(1, nrow(x$uni)), pch = if(is.numeric(x)) NULL else rep(19, nrow(x$uni)), ...)
{

   if(is.numeric(x)){
       y <- as.vector(x)
       names(y) <- names(x)  
       dotchart(y, xlab = axisLABEL, ...)
   }
   if(is.list(x)){
      if(length(col)==1) col <- rep(col, nrow(x$uni))
      if(length(pch)==1) pch <- rep(pch, nrow(x$uni))
      plot(x$q, x$uni[1, ], type = type, col = col[1], ylim = c(min(x$uni), max(x$uni)), pch = pch[1], , ylab=axisLABEL, xlab="q", ...)
      for(i in 1:nrow(x$uni)){
         lines(x$q, x$uni[i, ], type = type, col = col[i], pch = pch[i], ...)
      }
      if(legend[1]){
         legend(legendposi, legend = rownames(x$uni), col = col, lty = lty, pch = pch, ...) 
      }
   }    

}
