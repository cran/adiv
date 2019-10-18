plot.abgdivparam  <-
function(x, legend = TRUE, legendposi = "topright", type="b", col = if(is.numeric(x)) NULL else 1:nrow(x$div), lty = if(is.numeric(x)) NULL else rep(1, nrow(x$div)), pch = if(is.numeric(x)) NULL else 1:nrow(x$div), 
   ylim1 = range(x$div[c("Alpha","Gamma"), ]), ylim2 = NULL, ...)
{

   if(is.numeric(x)){
       y <- as.vector(x)
       names(y) <- names(x)  
       dotchart(y, xlab = "diversity", ...)
   }
   if(is.list(x)){
      if(nrow(x$div)==3){
      par(mfrow=c(2,1))
      plot(x$q, x$div[3, ], type = type, col = col[1], pch = pch[1], ylab="diversity", xlab="q", ylim = ylim1, ...)
      lines(x$q, x$div[1, ], type = type, col = col[2], pch = pch[2], ...)
      if(legend[1]){
         legend(legendposi, legend = rownames(x$div)[c(3,1)], col = col[1:2], lty = lty[1:2], pch = pch[1:2], ...) 
      }
      plot(x$q, x$div[2, ], type = type, col = col[3], pch = pch[3], ylab="diversity", xlab="q", ylim = ylim2,...)
      if(legend[1]){
         legend(legendposi, legend = rownames(x$div)[2], col = col[3], lty = lty[3], pch = pch[3], ...) 
      }
      }
      else{
      par(mfrow=c(3,1))
      plot(x$q, x$div[4, ], type = type, col = col[1], pch = pch[1], lty = lty[1], ylab="diversity", xlab="q", ylim = ylim1, ...)
      lines(x$q, x$div[1, ], type = type, col = col[2], pch = pch[2], lty = lty[2], ...)
      if(legend[1]){
         legend(legendposi, legend = rownames(x$div)[c(4,1)], col = col[1:2], lty = lty[1:2], pch = pch[1:2], ...) 
      }
      plot(x$q, x$div[2, ], type = type, col = col[3], pch = pch[3], lty = lty[3], ylab="diversity", xlab="q", ylim = ylim2,...)
      if(legend[1]){
         legend(legendposi, legend = rownames(x$div)[2], col = col[3], lty = lty[3], pch = pch[3], ...) 
      }    
      plot(x$q, x$div[3, ], type = type, col = col[4], pch = pch[4], lty = lty[4], ylab="dissimilarity", xlab="q", ylim = c(0, 1),...)
      if(legend[1]){
         legend(legendposi, legend = rownames(x$div)[3], col = col[4], lty = lty[4], pch = pch[4], ...) 
      }  
      }
   }
   par(mfrow=c(1,1))    
}
