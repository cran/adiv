plot.decdiv <- 
function (x, ncom = 1, col = "black", csize = 1, legend = TRUE, 
    ...) 
{
    if (!ncom %in% (1:ncol(x))) 
        stop("ncom must be the number of the column of x that must be plotted")
    if (!inherits(x, "decdiv")) 
        stop("x must be of class decdiv")
    phyape <- attributes(x)$phyl
    plot(phyape, plot = FALSE, ...)
    plotinfo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    MIN <- plotinfo$y.lim[1]
    MAX <- plotinfo$y.lim[2]   
    plot(phyape, y.lim = c(MIN-(MAX-MIN)/3.5, MAX), ...)
    plotinfo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    XX <- plotinfo$xx[-(1:length(phyape$tip.label))]
    YY <- plotinfo$yy[-(1:length(phyape$tip.label))]
    He <- min(c(max(YY) - min(YY), max(XX) - min(XX))) * csize
    he <- He/10
    sq <- x[phyape$node.label, ncom]/max(x[, ncom])
    radi <- sq * he
    symbols(XX, YY, circles = radi, bg = col, fg = "black", 
        inches = FALSE, add = TRUE)
    if (legend == TRUE) {
        parmar <- par()$mar
        par(mar = c(0.1, 0.1, 4.1, 2.1))
        sq0 <- pretty(c(x[phyape$node.label, ncom]), 4)
        l0 <- length(sq0)
        sq0 <- (sq0[1:(l0 - 1)] + sq0[2:l0])/2
        br0 <- (sq0/max(x[, ncom]) * he)
        sq0 <- round(sq0, digits = 3)
        cha <- as.character(sq0[1])
        for (i in (2:(length(sq0)))) cha <- paste(cha, sq0[i], 
            sep = " ")
        cex0 <- par("cex")
        yh <- max(c(strheight(cha, cex = cex0), br0))
        h <- strheight(cha, cex = cex0)
        y0 <- par("usr")[3] + yh/(par("usr")[2] - 
            par("usr")[1]) * (par("usr")[4] - par("usr")[3])
        x0 <- par("usr")[1] + h/2
        for (i in (1:(length(br0)))) {
            cha <- sq0[i]
            cha <- paste(" ", cha, sep = "")
            xh <- strwidth(cha, cex = cex0)
            text(x0 + xh/2, y0, cha, cex = cex0)
            z0 <- br0[i]
            x0 <- x0 + xh + z0
            symbols(x0, y0, circles = z0, bg = col, fg = "black", 
                add = TRUE, inches = FALSE)
            x0 <- x0 + z0
        }
        par(mar = parmar)
    }
}
