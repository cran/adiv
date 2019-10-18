plot.eveparam <-
function(x, legend = TRUE, legendposi = "topright", axisLABEL = "Evenness", type="b", col = if(is.numeric(x)) NULL else sample(colors(distinct = TRUE), nrow(x$eve)), lty = if(is.numeric(x)) NULL else rep(1, nrow(x$eve)), pch = if(is.numeric(x)) NULL else rep(19, nrow(x$eve)), ...)
{

   if(is.numeric(x)){
       y <- as.vector(x)
       names(y) <- names(x)  
       dotchart(y, xlab = axisLABEL, ...)
   }
   if(is.list(x)){
      if(length(col)==1) col <- rep(col, nrow(x$eve))
      if(length(pch)==1) pch <- rep(pch, nrow(x$eve))
      plot(x$q, x$eve[1, ], type = type, col = col[1], ylim = c(min(x$eve), max(x$eve)), pch = pch[1], , ylab = axisLABEL, xlab="q", ...)
      for(i in 1:nrow(x$eve)){
         lines(x$q, x$eve[i, ], type = type, col = col[i], pch = pch[i], ...)
      }
      if(legend[1]){
         legend(legendposi, legend = rownames(x$eve), col = col, lty = lty, pch = pch, ...) 
      }
   }    

}
