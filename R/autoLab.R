.autoLab <- function (x, y = NULL, labels = seq(along = x), cex = 1, method = c("SANN", 
    "GA"), allowSmallOverlap = FALSE, trace = FALSE, shadotext = FALSE, 
    doPlot = TRUE, ...) 
{
    if (!missing(y) && (is.character(y) || is.expression(y))) {
        labels <- y
        y <- NULL
    }
    if (is.factor(labels)) 
        labels <- as.character(labels)
    z <- xy.coords(x, y, recycle = TRUE)
    x <- z$x
    y <- z$y
    if (length(labels) < length(x)) 
        labels <- rep(labels, length(x))
    method <- match.arg(method)
    boundary <- par()$usr
    image_width <- boundary[2] - boundary[1]
    image_height <- boundary[4] - boundary[3]
    if (allowSmallOverlap) 
        nudgeFactor <- 0.02 * (abs(boundary[1] + (0+1i) * boundary[2] - 
            boundary[3] - (0+1i) * boundary[4]))
    n_labels <- length(x)
    xBoundary <- image_width * 0.01
    yBoundary <- image_height * 0.01
    width <- strwidth(labels, units = "user", cex = cex) + xBoundary
    height <- strheight(labels, units = "user", cex = cex) + 
        yBoundary
    gen_offset <- function(code) c(-1, -1, -1, 0, 0, 1, 1, 1)[code] * 
        (width/2) + (0+1i) * c(-1, 0, 1, -1, 1, -1, 0, 1)[code] * 
        (height/2)
    rect_intersect <- function(xy1, offset1, xy2, offset2) {
        w <- pmin(Re(xy1 + offset1/2), Re(xy2 + offset2/2)) - 
            pmax(Re(xy1 - offset1/2), Re(xy2 - offset2/2))
        h <- pmin(Im(xy1 + offset1/2), Im(xy2 + offset2/2)) - 
            pmax(Im(xy1 - offset1/2), Im(xy2 - offset2/2))
        w[w <= 0] <- 0
        h[h <= 0] <- 0
        w * h
    }
    nudge <- function(offset) {
        doesIntersect <- rect_intersect(xy[rectidx1] + offset[rectidx1], 
            rectv[rectidx1], xy[rectidx2] + offset[rectidx2], 
            rectv[rectidx2]) > 0
        pyth <- abs(xy[rectidx1] + offset[rectidx1] - xy[rectidx2] - 
            offset[rectidx2])/nudgeFactor
        eps <- 1e-10
        for (i in which(doesIntersect & pyth > eps)) {
            idx1 <- rectidx1[i]
            idx2 <- rectidx2[i]
            vect <- (xy[idx1] + offset[idx1] - xy[idx2] - offset[idx2])/pyth[idx1]
            offset[idx1] <- offset[idx1] + vect
            offset[idx2] <- offset[idx2] - vect
        }
        offset
    }
    objective <- function(gene) {
        offset <- gen_offset(gene)
        if (allowSmallOverlap) 
            offset <- nudge(offset)
        if (!is.null(rectidx1)) 
            area <- sum(rect_intersect(xy[rectidx1] + offset[rectidx1], 
                rectv[rectidx1], xy[rectidx2] + offset[rectidx2], 
                rectv[rectidx2]))
        else area <- 0
        n_outside <- sum(Re(xy + offset - rectv/2) < boundary[1] | 
            Re(xy + offset + rectv/2) > boundary[2] | Im(xy + 
            offset - rectv/2) < boundary[3] | Im(xy + offset + 
            rectv/2) > boundary[4])
        area + n_outside * image_width * image_height
    }
    xy <- x + (0+1i) * y
    rectv <- width + (0+1i) * height
    rectidx1 <- rectidx2 <- array(0, (length(x)^2 - length(x))/2)
    k <- 0
    for (i in 1:length(x)) for (j in seq(len = (i - 1))) {
        k <- k + 1
        rectidx1[k] <- i
        rectidx2[k] <- j
    }
    canIntersect <- rect_intersect(xy[rectidx1], 2 * rectv[rectidx1], 
        xy[rectidx2], 2 * rectv[rectidx2]) > 0
    rectidx1 <- rectidx1[canIntersect]
    rectidx2 <- rectidx2[canIntersect]
    if (trace) 
        cat("possible intersects =", length(rectidx1), "\n")
    if (trace) 
        cat("portion covered =", sum(rect_intersect(xy, rectv, 
            xy, rectv))/(image_width * image_height), "\n")
    SANN <- function() {
        gene <- rep(1, n_labels)
        gene <- gene + as.integer(Im(xy) > min(Im(xy)) + diff(range(Im(xy)))/3)
        gene <- gene + as.integer(Im(xy) > min(Im(xy)) + 2 * 
            diff(range(Im(xy)))/3)
        gene <- gene + 3 * as.integer(Re(xy) > min(Re(xy)) + 
            diff(range(Re(xy)))/3)
        gene <- gene + 3 * as.integer(Re(xy) > min(Re(xy)) + 
            2 * diff(range(Re(xy)))/3)
        gene[gene >= 6] <- gene[gene >= 6] - 1
        score <- objective(gene)
        bestgene <- gene
        bestscore <- score
        T <- 2.5
        for (i in 1:50) {
            k <- 1
            for (j in 1:50) {
                newgene <- gene
                newgene[sample(1:n_labels, 1)] <- sample(1:8, 
                  1)
                newscore <- objective(newgene)
                if (newscore < score || runif(1) < 1 - exp((newscore - 
                  score)/T)) {
                  k <- k + 1
                  score <- newscore
                  gene <- newgene
                }
                if (score <= bestscore) {
                  bestscore <- score
                  bestgene <- gene
                }
                if (bestscore == 0 || k == 10) 
                  break
            }
            if (bestscore == 0) 
                break
            if (trace) 
                cat("overlap area =", bestscore, "\n")
            T <- 0.9 * T
        }
        if (trace) 
            cat("overlap area =", bestscore, "\n")
        nx <- Re(xy + gen_offset(bestgene))
        ny <- Im(xy + gen_offset(bestgene))
        list(x = nx, y = ny)
    }
    xy <- SANN()
    shadowtext <- function(xy, labels, col = "black", bg = "white", 
        theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, ...) {
        xy <- xy.coords(xy)
        xo <- r * strwidth("A")
        yo <- r * strheight("A")
        for (i in theta) text(xy$x + cos(i) * xo, xy$y + sin(i) * 
            yo, labels, col = bg, ...)
        text(xy$x, xy$y, labels, col = col, ...)
    }
    if (doPlot) 
        if (shadotext == TRUE) 
            shadowtext(xy, labels, cex = cex, ...)
        else text(xy, labels, cex = cex, ...)
    invisible(xy)
}
