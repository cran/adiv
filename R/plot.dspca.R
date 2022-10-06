plot.dspca <- function(x, xaxis = 1, yaxis = 2, labels = TRUE, arrows = TRUE, points = FALSE, autolab = TRUE, title = NULL, colors = NULL, type = c("X&Y","X","Y"), zoom = TRUE, ...) {

    if(!inherits(x, "dspca")) 
        stop("x must be of class dspca")
    labels <- labels[1]
    type <- type[1]
    if(!type%in%c("X&Y","X","Y")) 
        type="X&Y"
    nSp <- nrow(x$X)
    nCo <- nrow(x$Y)
    if(type == "X&Y") {
        par(mfrow = c(1,2))
        if(!is.null(title) & length(title)>1) {
            titreSp <- title[1]
            titreCo <- title[2]
        }
        else {
            titreSp <- "Species map"
            titreCo <- "Community map"
        }
        if(!is.null(colors) & is.list(colors) & length(colors)==2) {
            colSp <- colors[[1]]
            colCo <- colors[[2]]
        }
        else {
            colSp <- "black"
            colCo <- "black"
        }
        if(zoom) {
            xlimSp <- xlimCo <- range(c(x$X[, xaxis], x$Y[, xaxis]))*110/100
            ylimSp <- ylimCo <- range(c(x$X[, yaxis], x$Y[, yaxis]))*110/100
        }
        else {
            xlimSp <- c(-1.1,1.1)
            ylimSp <- c(-1.1,1.1)
            xlimCo <- c(-1.1,1.1)
            ylimCo <- c(-1.1,1.1)

        }
        plot(0, 0, main = titreSp, xlab = paste("PC", xaxis, sep=""), 
            ylab = paste("PC", yaxis, sep=""), xlim = xlimSp, ylim = ylimSp, asp = 1, ...)
        abline(h=0, lty = 2)
        abline(v=0, lty = 2)
        if(arrows)
            arrows(rep(0, nSp), rep(0,nSp), x$X[, xaxis], x$X[, yaxis], col = colSp, ...)
        if(points)
            points(x$X[, xaxis], x$X[, yaxis], col = colSp, ...)
        if(labels) {
            if(autolab) 
                 .autoLab(x$X[, xaxis], x$X[, yaxis], labels = rownames(x$X), col = colSp, ...)
            else {
               pos <- rep(2, nSp)
               pos[x$X[, xaxis] > 0] <- 4
               text(x$X[, xaxis], x$X[, yaxis], labels = rownames(x$X), col = colSp, pos = pos, ...)
            }
        }
        symbols(0,0,1, inches=F, add=TRUE)
        plot(0, 0, main = titreCo, xlab = paste("PC", xaxis, sep=""), 
            ylab = paste("PC", yaxis, sep=""), xlim = xlimCo, ylim = ylimCo, asp = 1, ...)
        abline(h=0, lty = 2)
        abline(v=0, lty = 2)
        if(arrows)
            arrows(rep(0, nCo), rep(0,nCo), x$Y[, xaxis], x$Y[, yaxis], col = colCo, ...)
        if(points)
            points(x$Y[, xaxis], x$Y[, yaxis], col = colCo, ...)
        if(labels) {
            if(autolab) 
                 .autoLab(x$Y[, xaxis], x$Y[, yaxis], labels = rownames(x$Y), col = colCo, ...)
            else {
               pos <- rep(2, nCo)
               pos[x$Y[, xaxis] > 0] <- 4
               text(x$Y[, xaxis], x$Y[, yaxis], labels = rownames(x$Y), col = colCo, pos = pos, ...)
            }
        }
        symbols(0,0,1, inches=F, add=TRUE)       
        par(mfrow=c(1,1))
    }
    if(type == "X") {
        if(!is.null(title))
            titreSp <- title[1]
        else
            titreSp <- "Species map"
        if(!is.null(colors) & is.list(colors))
            colSp <- colors[[1]]
        else
            colSp <- "black"
        if(zoom) {
            xlimSp <- range(x$X[, xaxis])*110/100
            ylimSp <- range(x$X[, yaxis])*110/100
        }
        else {
            xlimSp <- c(-1.1,1.1)
            ylimSp <- c(-1.1,1.1)
        }
        plot(0, 0, main = titreSp, xlab = paste("PC", xaxis, sep=""), 
            ylab = paste("PC", yaxis, sep=""), xlim = xlimSp, ylim = ylimSp, asp = 1, ...)
        abline(h=0, lty = 2)
        abline(v=0, lty = 2)
        if(arrows)
            arrows(rep(0, nSp), rep(0,nSp), x$X[, xaxis], x$X[, yaxis], col = colSp, ...)
        if(points)
            points(x$X[, xaxis], x$X[, yaxis], col = colSp, ...)
        if(labels) {
            if(autolab) 
                 .autoLab(x$X[, xaxis], x$X[, yaxis], labels = rownames(x$X), col = colSp, ...)
            else {
               pos <- rep(2, nSp)
               pos[x$X[, xaxis] > 0] <- 4
               text(x$X[, xaxis], x$X[, yaxis], labels = rownames(x$X), col = colSp, pos = pos, ...)
            }
        }
        symbols(0,0,1, inches=F, add=TRUE)
    }
    if(type == "Y") {
        if(!is.null(title) & length(title)>1)
            titreCo <- title[2]
        else if(!is.null(title) & length(title)==1)
            titreCo <- title[1]
        else
            titreCo <- "Community map"
        if(!is.null(colors) & is.list(colors) & length(colors)==2) 
            colCo <- colors[[2]]
        else if(!is.null(colors) & !is.list(colors) & is.vector(colors))
            colCo <- colors
        else
            colCo <- "black"
        if(zoom) {
            xlimCo <- range(x$Y[, xaxis])*110/100  
            ylimCo <- range(x$Y[, yaxis])*110/100
        }
        else {
            xlimCo <- c(-1.1,1.1)
            ylimCo <- c(-1.1,1.1)
        }
        plot(0, 0, main = titreCo, xlab = paste("PC", xaxis, sep=""), 
            ylab = paste("PC", yaxis, sep=""), xlim = xlimCo, ylim = ylimCo, asp = 1, ...)
        abline(h=0, lty = 2)
        abline(v=0, lty = 2)
        if(arrows)
            arrows(rep(0, nCo), rep(0,nCo), x$Y[, xaxis], x$Y[, yaxis], col = colCo, ...)
        if(points)
            points(x$Y[, xaxis], x$Y[, yaxis], col = colCo, ...)
        if(labels) {
            if(autolab) 
                 .autoLab(x$Y[, xaxis], x$Y[, yaxis], labels = rownames(x$Y), col = colCo, ...)
            else {
               pos <- rep(2, nCo)
               pos[x$Y[, xaxis] > 0] <- 4
               text(x$Y[, xaxis], x$Y[, yaxis], labels = rownames(x$Y), col = colCo, pos = pos, ...)
            }
        }
        symbols(0,0,1, inches=F, add=TRUE)       
        par(mfrow=c(1,1))
    }


}

