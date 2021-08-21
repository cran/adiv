plot.FPdivparam <-
function(x, legend = TRUE, legendposi = "topright", axisLABEL = "FP-diversity", type="b", col = if(is.numeric(x)) NULL else sample(colors(distinct = TRUE), nrow(x$div)), lty = if(is.numeric(x)) NULL else rep(1, nrow(x$div)), pch = if(is.numeric(x)) NULL else rep(19, nrow(x$div)), ...)
{

   if(is.data.frame(x)){
       y <- as.vector(x[,1])
       names(y) <- rownames(x)  
       dotchart(y, xlab = axisLABEL, ...)
   }
   else{
      if(length(col)==1) col <- rep(col, nrow(x$div))
      if(length(pch)==1) pch <- rep(pch, nrow(x$div))
      plot(x$palpha, x$div[1, ], type = type, col = col[1], ylim = c(min(x$div), max(x$div)), pch = pch[1], , ylab = axisLABEL, xlab="alpha", ...)
      for(i in 1:nrow(x$div)){
         lines(x$palpha, x$div[i, ], type = type, col = col[i], pch = pch[i], ...)
      }
      if(legend[1]){
         legend(legendposi, legend = rownames(x$div), col = col, lty = lty, pch = pch, ...) 
      }
   }    

}
