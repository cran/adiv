plot.phylo4d <- function (x, trait = names(tdata(p4d)), center = TRUE, scale = TRUE, 
    plot.type = "barplot", tree.ladderize = FALSE, tree.type = "phylogram", 
    tree.ratio = NULL, tree.xlim = NULL, tree.open.angle = 0, 
    tree.open.crown = TRUE, show.tip = TRUE, tip.labels = NULL, 
    tip.col = "black", tip.cex = 1, tip.font = 3, tip.adj = 0, 
    data.xlim = NULL, bar.lwd = 10, bar.col = "grey35", 
    show.data.axis = TRUE, dot.col = "black", dot.pch = 20, 
    dot.cex = 2, cell.col = topo.colors(100), show.color.scale = TRUE, 
    show.trait = TRUE, trait.labels = NULL, trait.col = "black", 
    trait.cex = 1, trait.font = 1, trait.bg.col = "grey90", 
    error.bar.sup = NULL, error.bar.inf = NULL, error.bar.col = 1, 
    show.box = FALSE, grid.vertical = TRUE, grid.horizontal = FALSE, 
    grid.col = "grey25", grid.lty = "dashed", ...) 
{

    p4d <- x

    orderGrArg <- function (x, n.tips, n.traits, new.order, tips, default) 
    {
    x.dep <- deparse(substitute(x))
    if (is.vector(x)) {
        if (is.null(names(x))) {
            x <- rep(x, length.out = n.tips)
            x <- x[new.order]
        }
        else {
            if (any(tips %in% names(x))) {
                y <- rep(default, n.tips)
                names(y) <- tips
                y[names(x)] <- x[names(x)]
                x <- y[tips]
            }
            else {
                stop(paste("Phylogenetic tips do not match with names of", 
                  x.dep))
            }
        }
    }
    if (is.matrix(x)) {
        if (is.null(rownames(x))) {
            x <- x[new.order, ]
        }
        else {
            if (any(tips %in% rownames(x))) {
                y <- matrix(default, nrow = nrow(x), ncol = ncol(x))
                rownames(y) <- tips
                y[rownames(x), ] <- x[rownames(x), ]
                x <- y[tips, ]
            }
            else {
                stop(paste("Phylogenetic tips do not match with row names of", 
                  x.dep))
            }
        }
    }
    x <- matrix(x, nrow = n.tips, ncol = n.traits)
    return(x)
    }

    layouterize <- function (n.traits, show.tip) 
    {
        if (show.tip) {
            res <- matrix(c(n.traits + 2, 1:(n.traits + 1)), nrow = 1)
        }
        else {
            res <- matrix(c(n.traits + 1, 1:(n.traits)), nrow = 1)
        }
        return(res)
    }

    layouterizeRatio <- function (tree.ratio, n.traits, show.tip) 
    {
        if (!is.null(tree.ratio)) {
            if (show.tip) {
                res <- c(tree.ratio, rep((1 - tree.ratio)/(n.traits + 
                1), n.traits + 1))
            }
            else {
                res <- c(tree.ratio, rep((1 - tree.ratio)/n.traits, 
                n.traits))
            }
        }
        else {
            if (show.tip) {
                res <- rep(1, n.traits + 1)
            }
            else {
                res <- rep(1, n.traits)
            }
        }
        return(res)
    }

    matchTipsAndTraits <- function (x, p4d = NULL, p4d.tips = NULL, p4d.traits = NULL, 
    subset = TRUE) 
    {
        if (!is.null(p4d) & is(p4d, "phylo4d")) {
            p4d.tips <- tipLabels(p4d)
            p4d.traits <- colnames(tdata(p4d))
        }
        if (!all(p4d.tips %in% rownames(x))) {
            stop("Rows names of !!! do not match with tips names")
        }
        if (!all(p4d.traits %in% colnames(x))) {
            stop("Columns names of !!! do not match with traits names")
        }
        if (subset) {
            x <- x[p4d.tips, p4d.traits]
        }
        return(x)
    }

    plotPhyloDisabled <- function (x, type = "phylogram", use.edge.length = TRUE, 
        node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE, 
        edge.color = "black", edge.width = 1, edge.lty = 1, 
        font = 3, cex = par("cex"), adj = NULL, srt = 0, no.margin = FALSE, 
        root.edge = FALSE, label.offset = 0, underscore = FALSE, 
        x.lim = NULL, y.lim = NULL, direction = "rightwards", 
        lab4ut = NULL, tip.color = "black", plot = TRUE, rotate.tree = 0, 
        open.angle = 0, node.depth = 1, ...) 
    {
        Ntip <- length(x$tip.label)
        if (Ntip < 2) {
            warning("found less than 2 tips in the tree")
            return(NULL)
        }
        if (any(tabulate(x$edge[, 1]) == 1)) 
            stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
        Nedge <- dim(x$edge)[1]
        Nnode <- x$Nnode
        if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
            stop("tree badly conformed; cannot plot. Check the edge matrix.")
        ROOT <- Ntip + 1
        type <- match.arg(type, c("phylogram", "cladogram", 
            "fan", "unrooted", "radial"))
        direction <- match.arg(direction, c("rightwards", "leftwards", 
            "upwards", "downwards"))
        if (is.null(x$edge.length)) 
            use.edge.length <- FALSE
        if (type %in% c("unrooted", "radial") || !use.edge.length || 
            is.null(x$root.edge) || !x$root.edge) 
            root.edge <- FALSE
        if (type == "fan" && root.edge) {
            warning("drawing root edge with type = 'fan' is not yet supported")
            root.edge <- FALSE
        }
        phyloORclado <- type %in% c("phylogram", "cladogram")
        horizontal <- direction %in% c("rightwards", "leftwards")
        xe <- x$edge
        if (phyloORclado) {
            phyOrder <- attr(x, "order")
            if (is.null(phyOrder) || phyOrder != "cladewise") {
                x <- reorder(x)
                if (!identical(x$edge, xe)) {
                    ereorder <- match(x$edge[, 2], xe[, 2])
                    if (length(edge.color) > 1) {
                      edge.color <- rep(edge.color, length.out = Nedge)
                      edge.color <- edge.color[ereorder]
                    }
                    if (length(edge.width) > 1) {
                      edge.width <- rep(edge.width, length.out = Nedge)
                      edge.width <- edge.width[ereorder]
                    }
                    if (length(edge.lty) > 1) {
                      edge.lty <- rep(edge.lty, length.out = Nedge)
                      edge.lty <- edge.lty[ereorder]
                    }
                }
            }
            yy <- numeric(Ntip + Nnode)
            TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
            yy[TIPS] <- 1:Ntip
        }
        z <- reorder(x, order = "postorder")
        if (phyloORclado) {
        if (is.null(node.pos)) 
            node.pos <- if (type == "cladogram" && !use.edge.length) 
                2
            else 1
        if (node.pos == 1) 
            yy <- node.height(x, clado.style = FALSE)
        else {
            xx <- node.depth(x) - 1
            yy <- node.height(x, clado.style = TRUE)
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- node.depth(x) - 1
                xx <- max(xx) - xx
            }
            else {
                xx <- node.depth.edgelength(x)
            }
        }
        else {
            twopi <- 2 * pi
            rotate.tree <- twopi * rotate.tree/360
            if (type != "unrooted") {
                TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
                xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                    length.out = Ntip)
                theta <- double(Ntip)
                theta[TIPS] <- xx
                theta <- c(theta, numeric(Nnode))
            }
            switch(type, fan = {
                theta <- node.height(x)
                if (use.edge.length) {
                    r <- node.depth.edgelength(x)
                } else {
                   r <- node.depth(x)
                    r <- 1/r
                }
                theta <- theta + rotate.tree
                xx <- r * cos(theta)
                yy <- r * sin(theta)
            }, unrooted = {
                nb.sp <- node.depth(x)
                XY <- if (use.edge.length) ape::unrooted.xy(Ntip, 
                    Nnode, z$edge, z$edge.length, nb.sp, rotate.tree) else ape::unrooted.xy(Ntip, 
                    Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
                xx <- XY$M[, 1] - min(XY$M[, 1])
                yy <- XY$M[, 2] - min(XY$M[, 2])
            }, radial = {
                r <- node.depth(x)
                r[r == 1] <- 0
                r <- 1 - r/Ntip
                theta <- node.height(x) + rotate.tree
                xx <- r * cos(theta)
                yy <- r * sin(theta)
            })
        }
        if (phyloORclado) {
            if (!horizontal) {
                tmp <- yy
                yy <- xx
                xx <- tmp - min(tmp) + 1
            }
            if (root.edge) {
                if (direction == "rightwards") 
                    xx <- xx + x$root.edge
                if (direction == "upwards") 
                    yy <- yy + x$root.edge
            }
        }
        if (no.margin) 
            par(mai = rep(0, 4))
        if (show.tip.label) 
            nchar.tip.label <- nchar(x$tip.label)
        max.yy <- max(yy)
        if (is.null(x.lim)) {
            if (phyloORclado) { 
               if (horizontal) {
                    x.lim <- c(0, NA)
                    pin1 <- par("pin")[1]
                    strWi <- strwidth(x$tip.label, "inches", 
                      cex = cex)
                    xx.tips <- xx[1:Ntip] * 1.04
                    alp <- try(uniroot(function(a) max(a * xx.tips + 
                      strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                      tmp <- max(xx.tips) * 1.5
                    else {
                      tmp <- if (show.tip.label) 
                        max(xx.tips + strWi/alp)
                      else max(xx.tips)
                    }
                    if (show.tip.label) 
                      tmp <- tmp + label.offset
                    x.lim[2] <- tmp
                }
                else x.lim <- c(1, Ntip)
            }
            else switch(type, fan = {
                if (show.tip.label) {
                    offset <- max(nchar.tip.label * 0.018 * max.yy * 
                      cex)
                    x.lim <- range(xx) + c(-offset, offset)
                } else x.lim <- range(xx)
            }, unrooted = {
                if (show.tip.label) {
                    offset <- max(nchar.tip.label * 0.018 * max.yy * 
                      cex)
                    x.lim <- c(0 - offset, max(xx) + offset)
                } else x.lim <- c(0, max(xx))
            }, radial = {
                if (show.tip.label) {
                    offset <- max(nchar.tip.label * 0.03 * cex)
                    x.lim <- c(-1 - offset, 1 + offset)
                } else x.lim <- c(-1, 1)
            })
        }
        else if (length(x.lim) == 1) {
            x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
                x.lim[1] <- 1
            if (type %in% c("fan", "unrooted") && show.tip.label) 
                x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                    cex)
            if (type == "radial") 
                x.lim[1] <- if (show.tip.label) 
                    -1 - max(nchar.tip.label * 0.03 * cex)
                else -1
        }
        if (phyloORclado && direction == "leftwards") 
            xx <- x.lim[2] - xx
        if (is.null(y.lim)) {
            if (phyloORclado) {
                if (horizontal) 
                    y.lim <- c(1, Ntip)
                else {
                    y.lim <- c(0, NA)
                    pin2 <- par("pin")[2]
                    strWi <- strwidth(x$tip.label, "inches", 
                      cex = cex)
                    yy.tips <- yy[1:Ntip] * 1.04
                    alp <- try(uniroot(function(a) max(a * yy.tips + 
                      strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
                    if (is.character(alp)) 
                      tmp <- max(yy.tips) * 1.5
                    else {
                      tmp <- if (show.tip.label) 
                        max(yy.tips + strWi/alp)
                      else max(yy.tips)
                    }
                    if (show.tip.label) 
                      tmp <- tmp + label.offset
                   y.lim[2] <- tmp
                }
            }
            else switch(type, fan = {
                if (show.tip.label) {
                    offset <- max(nchar.tip.label * 0.018 * max.yy * 
                  cex)
                    y.lim <- c(min(yy) - offset, max.yy + offset)
                } else y.lim <- c(min(yy), max.yy)
            }, unrooted = {
                if (show.tip.label) {
                    offset <- max(nchar.tip.label * 0.018 * max.yy * 
                      cex)
                    y.lim <- c(0 - offset, max.yy + offset)
                } else y.lim <- c(0, max.yy)
            }, radial = {
                if (show.tip.label) {
                    offset <- max(nchar.tip.label * 0.03 * cex)
                    y.lim <- c(-1 - offset, 1 + offset)
                } else y.lim <- c(-1, 1)
            })
        }
        else if (length(y.lim) == 1) {
            y.lim <- c(0, y.lim)
            if (phyloORclado && horizontal) 
                y.lim[1] <- 1
            if (type %in% c("fan", "unrooted") && show.tip.label) 
                y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                    cex)
            if (type == "radial") 
                y.lim[1] <- if (show.tip.label) 
                    -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
                else -1
        }
        if (phyloORclado && direction == "downwards") 
            yy <- y.lim[2] - yy
        if (phyloORclado && root.edge) {
            if (direction == "leftwards") 
                x.lim[2] <- x.lim[2] + x$root.edge
            if (direction == "downwards") 
                y.lim[2] <- y.lim[2] + x$root.edge
        }
        asp <- if (type %in% c("fan", "radial", "unrooted")) 
            1
        else NA
        L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
        assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, 
            yy = yy)), envir = ape::.PlotPhyloEnv)
        return(L)
    }

    p4 <- phylobase::extractTree(p4d)
    phy <- as(p4, "phylo")
    if (tree.ladderize) {
        phy <- ladderize(phy)
    }
    new.order <- phy$edge[, 2][!phy$edge[, 2] %in% phy$edge[, 
        1]]
    tips <- phy$tip.label[new.order]
    n.tips <- length(tips)
    X <- tdata(p4d, type = "tip")
    X <- X[new.order, trait]
    X <- scale(X, center = center, scale = scale)
    X <- as.data.frame(X)
    colnames(X) <- trait
    n.traits <- ncol(X)
    if (is.numeric(trait)) {
        trait <- names(tdata(p4d))[trait]
    }
    tree.type <- match.arg(tree.type, c("phylogram", "cladogram", 
        "fan"))
    plot.type <- match.arg(plot.type, c("barplot", "dotplot", 
        "gridplot"))
    if (!is.null(error.bar.inf)) {
        error.bar.inf <- as.matrix(error.bar.inf)
        error.bar.inf[is.na(error.bar.inf)] <- 0
        error.bar.inf <- matchTipsAndTraits(error.bar.inf, p4d.tips = tips, 
            p4d.traits = trait)
    }
    if (!is.null(error.bar.sup)) {
        error.bar.sup <- as.matrix(error.bar.sup)
        error.bar.sup[is.na(error.bar.sup)] <- 0
        error.bar.sup <- matchTipsAndTraits(error.bar.sup, p4d.tips = tips, 
            p4d.traits = trait)
    }
    if (!is.null(error.bar.inf)) {
        arrow.inf <- X
        arrow.inf[X > 0] <- X[X > 0] - error.bar.inf[X > 0]
        arrow.inf[X < 0] <- X[X < 0] + error.bar.inf[X < 0]
    }
    if (!is.null(error.bar.sup)) {
        arrow.sup <- X
        arrow.sup[X > 0] <- X[X > 0] + error.bar.sup[X > 0]
        arrow.sup[X < 0] <- X[X < 0] - error.bar.sup[X < 0]
    }
    if (is.null(data.xlim)) {
        if (center & scale) {
            data.xlim <- matrix(rep(NA, n.traits * 2), nrow = 2, 
                dimnames = list(c("xlim.min", "xlim.max"), 
                  trait))
            if (!is.null(error.bar.inf) & !is.null(error.bar.sup)) {
                data.xlim[1, ] <- floor(min(arrow.inf, arrow.sup, 
                  na.rm = TRUE))
                data.xlim[2, ] <- ceiling(max(arrow.inf, arrow.sup, 
                  na.rm = TRUE))
            }
            else if (!is.null(error.bar.inf)) {
                data.xlim[1, ] <- floor(min(arrow.inf, na.rm = TRUE))
                data.xlim[2, ] <- ceiling(max(arrow.inf, na.rm = TRUE))
            }
            else if (!is.null(error.bar.sup)) {
                data.xlim[1, ] <- floor(min(arrow.sup, na.rm = TRUE))
                data.xlim[2, ] <- ceiling(max(arrow.sup, na.rm = TRUE))
            }
            else {
                data.xlim[1, ] <- floor(min(X, na.rm = TRUE))
                data.xlim[2, ] <- ceiling(max(X, na.rm = TRUE))
            }
        }
        else {
            data.xlim <- matrix(NA, nrow = 2, ncol = n.traits, 
                dimnames = list(c("xlim.min", "xlim.max"), 
                  trait))
            data.xlim[1, ] <- apply(X, 2, min, na.rm = TRUE)
            data.xlim[1, apply(X, 2, min, na.rm = TRUE) * apply(X, 
                2, max, na.rm = TRUE) > 0 & apply(X, 2, min, 
                na.rm = TRUE) > 0] <- 0
            data.xlim[2, ] <- apply(X, 2, max, na.rm = TRUE)
            data.xlim[2, apply(X, 2, min, na.rm = TRUE) * apply(X, 
                2, max, na.rm = TRUE) > 0 & apply(X, 2, max, 
                na.rm = TRUE) < 0] <- 0
            if (!is.null(error.bar.inf) & !is.null(error.bar.sup)) {
                data.xlim[1, ] <- apply(cbind(apply(arrow.inf, 
                  2, min, na.rm = TRUE), apply(arrow.sup, 2, 
                  min)), 1, min, na.rm = TRUE)
                data.xlim[2, ] <- apply(cbind(apply(arrow.inf, 
                  2, max, na.rm = TRUE), apply(arrow.sup, 2, 
                  max)), 1, max, na.rm = TRUE)
            }
            else {
                if (!is.null(error.bar.inf)) {
                  data.xlim[1, ] <- apply(cbind(apply(arrow.inf, 
                    2, min, na.rm = TRUE), data.xlim[1, ]), 1, 
                    min, na.rm = TRUE)
                  data.xlim[2, ] <- apply(cbind(apply(arrow.inf, 
                    2, max, na.rm = TRUE), data.xlim[2, ]), 1, 
                    max, na.rm = TRUE)
                }
                if (!is.null(error.bar.sup)) {
                  data.xlim[1, ] <- apply(cbind(apply(arrow.sup, 
                    2, min, na.rm = TRUE), data.xlim[1, ]), 1, 
                    min, na.rm = TRUE)
                  data.xlim[2, ] <- apply(cbind(apply(arrow.sup, 
                    2, max, na.rm = TRUE), data.xlim[2, ]), 1, 
                    max, na.rm = TRUE)
                }
            }
        }
    }
    else if (is.vector(data.xlim) & length(data.xlim) == 2) {
        data.xlim <- matrix(rep(data.xlim, n.traits), nrow = 2, 
            dimnames = list(c("xlim.min", "xlim.max"), 
                trait))
    }
    else if (is.matrix(data.xlim)) {
        if (isTRUE(all.equal(dim(data.xlim), c(2, n.traits)))) {
            rownames(data.xlim) <- c("xlim.min", "xlim.max")
            colnames(data.xlim) <- trait
        }
        else {
            stop("Invalid 'data.xlim' argument: wrong matrix dimensions")
        }
    }
    else {
        stop("Invalid 'data.xlim' argument")
    }
    ylim <- c(1, n.tips)
    if (plot.type == "barplot") {
        bar.col <- orderGrArg(bar.col, n.tips = n.tips, n.traits = n.traits, 
            new.order = new.order, tips = tips, default = "grey35")
    }
    if (plot.type == "dotplot") {
        dot.col <- orderGrArg(dot.col, n.tips = n.tips, n.traits = n.traits, 
            new.order = new.order, tips = tips, default = 1)
        dot.pch <- orderGrArg(dot.pch, n.tips = n.tips, n.traits = n.traits, 
            new.order = new.order, tips = tips, default = 1)
        dot.cex <- orderGrArg(dot.cex, n.tips = n.tips, n.traits = n.traits, 
            new.order = new.order, tips = tips, default = 1)
    }
    if (!is.null(error.bar.inf) | !is.null(error.bar.sup)) {
        error.bar.col <- orderGrArg(error.bar.col, n.tips = n.tips, 
            n.traits = n.traits, new.order = new.order, tips = tips, 
            default = 1)
    }
    if (is.null(tip.labels)) {
        tip.labels <- tips
    }
    else {
        tip.labels <- orderGrArg(tip.labels, n.tips = n.tips, 
            n.traits = n.traits, new.order = new.order, tips = tips, 
            default = "")
    }
    tip.col <- orderGrArg(tip.col, n.tips = n.tips, n.traits = n.traits, 
        new.order = new.order, tips = tips, default = 1)
    tip.cex <- orderGrArg(tip.cex, n.tips = n.tips, n.traits = n.traits, 
        new.order = new.order, tips = tips, default = 1)
    tip.font <- orderGrArg(tip.font, n.tips = n.tips, n.traits = n.traits, 
        new.order = new.order, tips = tips, default = 3)
    if (is.null(trait.labels)) {
        trait.labels <- trait
    }
    trait.labels <- rep(trait.labels, length.out = n.traits)
    trait.col <- rep(trait.col, length.out = n.traits)
    trait.cex <- rep(trait.cex, length.out = n.traits)
    trait.font <- rep(trait.font, length.out = n.traits)
    trait.bg.col <- rep(trait.bg.col, length.out = n.traits)
    if (is.null(tree.xlim)) {
        tree.xlim <- plotPhyloDisabled(phy, type = tree.type, 
            show.tip.label = FALSE, x.lim = NULL, y.lim = NULL, 
            no.margin = FALSE, direction = "rightwards", 
            plot = FALSE, ...)$x.lim
    }
    par.mar0 <- par("mar")
    par.lend0 <- par("lend")
    par.xpd0 <- par("xpd")
    if (tree.type == "phylogram" | tree.type == "cladogram") {
        if (plot.type %in% c("barplot", "dotplot")) {
            lay <- layouterize(n.traits = n.traits, show.tip = show.tip)
            lay.w <- layouterizeRatio(tree.ratio = tree.ratio, 
                n.traits = n.traits, show.tip = show.tip)
        }
        if (plot.type == "gridplot") {
            lay <- layouterize(n.traits = 1, show.tip = show.tip)
            lay.w <- layouterizeRatio(tree.ratio = tree.ratio, 
                n.traits = 1, show.tip = show.tip)
        }
        layout(lay, widths = lay.w)
        par(xpd = FALSE, mar = c(5, 1, 4, 0), lend = 1)
        fig.traits <- vector("list", n.traits)
        names(fig.traits) <- trait
        if (plot.type %in% c("barplot", "dotplot")) {
            for (i in 1:n.traits) {
                plot.new()
                plot.window(xlim = data.xlim[, i], ylim = ylim)
                fig.traits[[i]] <- par("fig")
                rect(par("usr")[1], par("usr")[3] - 
                  (3 * par("cxy")[2]), par("usr")[2], 
                  par("usr")[4], col = trait.bg.col[i], 
                  border = NA, xpd = TRUE)
                if (show.box) {
                  box()
                }
                if (grid.vertical) {
                  grid(NULL, NA, col = grid.col, lty = grid.lty)
                  abline(v = 0, lty = "solid", col = grid.col)
                }
                else {
                  abline(v = 0, lty = "solid", col = grid.col)
                }
                if (grid.horizontal) {
                  abline(h = 1:n.tips, col = grid.col, lty = grid.lty)
                }
                if (plot.type == "barplot") {
                  segments(x0 = 0, x1 = X[, i], y0 = 1:n.tips, 
                    lwd = bar.lwd, col = bar.col[, i])
                }
                if (plot.type == "dotplot") {
                  points(x = X[, i], y = 1:n.tips, col = dot.col[, 
                    i], pch = dot.pch[, i], cex = dot.cex[, i])
                }
                options(warn = -1)
                if (!is.null(error.bar.inf)) {
                  arrows(x0 = X[, i], x1 = arrow.inf[, i], y0 = 1:n.tips, 
                    lwd = 1, col = error.bar.col, angle = 90, 
                    length = 0.04)
                }
                if (!is.null(error.bar.sup)) {
                  arrows(x0 = X[, i], x1 = arrow.sup[, i], y0 = 1:n.tips, 
                    lwd = 1, col = error.bar.col, angle = 90, 
                    length = 0.04)
                }
                options(warn = 1)
                if (show.data.axis) {
                  axis(1)
                }
                if (show.trait) {
                  mtext(trait.labels[i], side = 1, line = 3, 
                    las = par("las"), col = trait.col[i], 
                    cex = trait.cex[i], font = trait.font[i])
                }
            }
        }
        if (plot.type == "gridplot") {
            plot.new()
            rect(par("usr")[1], par("usr")[3] - (3 * 
                par("cxy")[2]), par("usr")[2], par("usr")[4], 
                col = trait.bg.col[1], border = NA, xpd = TRUE)
            data.xlim[1, ] <- 0
            data.xlim[2, ] <- n.traits
            plot.window(xlim = data.xlim[, 1], ylim = ylim)
            fig.traits[[1]] <- par("fig")
            image(x = 0:n.traits, y = 1:n.tips, z = t(X), col = cell.col, 
                add = TRUE, xlab = "", ylab = "", 
                yaxs = FALSE, xaxs = FALSE)
            if (show.box) {
                box()
            }
            if (grid.horizontal) {
                abline(h = seq(1.5, n.tips - 0.5), col = grid.col, 
                  lty = grid.lty)
            }
            if (grid.vertical) {
                abline(v = seq(1, n.traits - 1), col = grid.col, 
                  lty = grid.lty)
            }
            if (show.trait) {
                mtext(trait.labels, at = seq(0.5, (n.traits - 
                  0.5)), side = 1, line = 1, las = par("las"), 
                  col = trait.col, cex = trait.cex, font = trait.font)
            }
        }
        if (show.tip) {
            plot.new()
            tip.xlim <- c(-1, 1)
            if (tip.adj < 0.5) 
                tip.xlim[1] <- -tip.adj/0.5
            if (tip.adj > 0.5) 
                tip.xlim[2] <- -2 * tip.adj + 2
            plot.window(xlim = tip.xlim, ylim = ylim)
            text(x = 0, y = 1:n.tips, labels = tip.labels, adj = tip.adj, 
                col = tip.col, cex = tip.cex, font = tip.font)
            fig.tip <- par("fig")
        }
        else {
            fig.tip <- NULL
            tip.xlim <- NULL
        }
        plot.phylo(phy, type = tree.type, show.tip.label = FALSE, 
            x.lim = tree.xlim, y.lim = NULL, no.margin = FALSE, 
            direction = "rightwards", ...)
        fig.tree <- par("fig")
        if (plot.type == "gridplot" & show.color.scale) {
            par(new = TRUE)
            plt.init <- par("plt")
            par(plt = c(par("plt")[1] + 0.05, par("plt")[2] - 
                0.2, 0.07, 0.1))
            plot.new()
            breaks <- seq(min(X), max(X), length.out = (length(cell.col) + 
                1))
            scale.xlim <- range(breaks)
            scale.ylim <- c(0, 1)
            plot.window(xlim = scale.xlim, ylim = scale.ylim)
            for (i in 1:length(cell.col)) {
                polygon(c(breaks[i], breaks[i + 1], breaks[i + 
                  1], breaks[i]), c(0, 0, 1, 1), col = cell.col[i], 
                  border = NA)
            }
            axis(1)
            par(plt = plt.init)
        }
        assign("last_barplotp4d", list(plot.type = plot.type, 
            show.tip = show.tip, layout = lay, fig.tree = fig.tree, 
            fig.traits = fig.traits, fig.tip = fig.tip, tree.xlim = tree.xlim, 
            data.xlim = data.xlim, tip.xlim = tip.xlim, ylim = ylim, 
            par.mar0 = par.mar0), envir = ape::.PlotPhyloEnv)
        layout(1)
    }
    if (tree.type == "fan") {
        par(lend = 1)
        if (is.null(tree.ratio)) {
            if (show.tip) {
                tree.ratio <- 1/(n.traits + 2)
            }
            else {
                tree.ratio <- 1/(n.traits + 1)
            }
        }
        plot.phylo(phy, type = tree.type, show.tip.label = FALSE, 
            x.lim = tree.xlim * (1/tree.ratio), y.lim = NULL, 
            no.margin = TRUE, open.angle = tree.open.angle, rotate.tree = 0, 
            ...)
        lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
        length.phylo <- max(sqrt(lp$xx^2 + lp$yy^2))
        if (show.tip) {
            length.gr0 <- (min(par("usr")[2] - par("usr")[1], 
                par("usr")[4] - par("usr")[3])/2 - 
                length.phylo)/(n.traits + 1)
        }
        else {
            length.gr0 <- (min(par("usr")[2] - par("usr")[1], 
                par("usr")[4] - par("usr")[3])/2 - 
                length.phylo)/n.traits
        }
        length.intergr <- 0.2 * length.gr0
        length.gr <- length.gr0 - length.intergr
        theta <- atan2(lp$xx[1:n.tips], lp$yy[1:n.tips])[new.order]
        theta[theta > (pi/2)] <- -pi - (pi - theta[theta > (pi/2)])
        cos.t <- cos(pi/2 - theta)
        sin.t <- sin(pi/2 - theta)
        theta.real.open <- diff(c(min(theta), max(theta)))
        real.open <- theta.real.open * 180/pi
        if (tree.open.crown) {
            theta.soft <- pi/2 - seq(-5, real.open + 5, length.out = 300) * 
                pi/180
        }
        else {
            theta.soft <- pi/2 - seq(0, 360) * pi/180 + 10 * 
                pi/180
        }
        cos.tsoft <- cos(pi/2 - theta.soft)
        sin.tsoft <- sin(pi/2 - theta.soft)
        for (i in 1:n.traits) {
            length.ring1 <- length.phylo + length.intergr * i + 
                length.gr * (i - 1) - 0.3 * length.intergr
            length.ring2 <- length.phylo + length.intergr * i + 
                length.gr * i + 0.3 * length.intergr
            xx1 <- length.ring1 * cos.tsoft
            xx2 <- length.ring2 * cos.tsoft
            yy1 <- length.ring1 * sin.tsoft
            yy2 <- length.ring2 * sin.tsoft
            polygon(c(xx1, rev(xx2)), c(yy1, rev(yy2)), col = trait.bg.col[i], 
                border = NA)
            if (abs(sign(min(data.xlim[, i])) + sign(max(data.xlim[, 
                i]))) == 2) {
                scaling.factor <- length.gr/max(abs(min(data.xlim[, 
                  i])), abs(max(data.xlim[, i])))
            }
            else {
                scaling.factor <- length.gr/diff(c(min(data.xlim[, 
                  i]), max(data.xlim[, i])))
            }
            X.scale <- X[, i] * scaling.factor
            data.xlim.scale <- data.xlim[, i] * scaling.factor
            if (!is.null(error.bar.inf)) {
                arrow.inf.scale <- arrow.inf[, i] * scaling.factor
            }
            if (!is.null(error.bar.sup)) {
                arrow.sup.scale <- arrow.sup[, i] * scaling.factor
            }
            if (plot.type == "barplot" | plot.type == "dotplot") {
                length.baseline <- (length.phylo + length.intergr * 
                  i + length.gr * (i - 1) + ifelse(min(data.xlim.scale) < 
                  0, abs(min(data.xlim.scale)), 0))
                length.baseline <- rep(length.baseline, n.tips)
                length.values <- length.baseline + X.scale
                if (!is.null(error.bar.inf)) {
                  length.arrow.inf <- length.baseline + arrow.inf.scale
                }
                if (!is.null(error.bar.sup)) {
                  length.arrow.sup <- length.baseline + arrow.sup.scale
                }
                length.baseline <- rep(length.baseline, length.out = length(cos.tsoft))
                lines(length.baseline * cos.tsoft, length.baseline * 
                  sin.tsoft, lwd = 1)
                if (grid.horizontal) {
                  segments(x0 = length.ring1 * cos.t, x1 = length.ring2 * 
                    cos.t, y0 = length.ring1 * sin.t, y1 = length.ring2 * 
                    sin.t, col = grid.col, lty = grid.lty)
                }
                if ((show.data.axis | grid.vertical)) {
                  if (tree.open.crown) {
                    theta.ax <- theta.soft[1]
                  }
                  else {
                    if (show.trait) {
                      theta.ax <- theta[1] + (360 - real.open) * 
                        (1/3) * pi/180
                    }
                    else {
                      theta.ax <- theta[1] + (360 - real.open) * 
                        (1/2) * pi/180
                    }
                  }
                  cos.tax <- cos(pi/2 - theta.ax)
                  sin.tax <- sin(pi/2 - theta.ax)
                  nint.ticks <- round((length.gr/min(par("usr")[2] - 
                    par("usr")[1], par("usr")[4] - 
                    par("usr")[3]))/3 * 100) - 1
                  if (min(data.xlim.scale) <= 0 & max(data.xlim.scale) >= 
                    0) {
                    ticks <- axisTicks(c(min(data.xlim.scale)/scaling.factor, 
                      max(data.xlim.scale)/scaling.factor), log = FALSE, 
                      nint = nint.ticks)
                  }
                  else {
                    if (abs(min(data.xlim.scale)) > max(data.xlim.scale)) {
                      ticks <- axisTicks(c(0, min(data.xlim.scale)/scaling.factor), 
                        log = FALSE, nint = nint.ticks)
                    }
                    else {
                      ticks <- axisTicks(c(0, max(data.xlim.scale)/scaling.factor), 
                        log = FALSE, nint = nint.ticks)
                    }
                  }
                  ticks <- ifelse(ticks > max(data.xlim.scale)/scaling.factor, 
                    NA, ticks)
                  length.ticks <- length.baseline[1] + ticks * 
                    scaling.factor
                  if (grid.vertical) {
                    for (j in 1:length(length.ticks)) {
                      lines(length.ticks[j] * cos.tsoft, length.ticks[j] * 
                        sin.tsoft, col = grid.col, lty = grid.lty)
                    }
                  }
                  if (show.data.axis) {
                    segments(x0 = length.ring1 * cos.tax, x1 = length.ring2 * 
                      cos.tax, y0 = length.ring1 * sin.tax, y1 = length.ring2 * 
                      sin.tax, lwd = 20, col = trait.bg.col[i])
                    text(x = length.ticks * cos.tax, y = length.ticks * 
                      sin.tax, labels = ticks, cex = tip.cex)
                  }
                }
            }
            if (plot.type == "barplot") {
                length.baseline <- rep(length.baseline, length.out = length(cos.t))
                segments(x0 = length.baseline * cos.t, x1 = length.values * 
                  cos.t, y0 = length.baseline * sin.t, y1 = length.values * 
                  sin.t, lwd = bar.lwd, col = bar.col[, i])
            }
            if (plot.type == "dotplot") {
                length.baseline <- rep(length.baseline, length.out = length(cos.t))
                points(x = length.values * cos.t, y = length.values * 
                  sin.t, col = dot.col[, i], pch = dot.pch[, 
                  i], cex = dot.cex[, i])
            }
            if (plot.type == "gridplot") {
                nc <- length(cell.col)
                X.cut <- as.numeric(cut(as.matrix(X), nc))
                grid.colors <- cell.col[X.cut]
                grid.colors <- matrix(grid.colors, ncol = n.traits)
                theta.grid1 <- theta + pi/2 - (theta[1] + theta[2])/2
                theta.grid2 <- theta - pi/2 + (theta[1] + theta[2])/2
                for (k in 1:length(theta)) {
                  tile.x1 <- length.ring1 * cos(pi/2 - seq(theta.grid1[k], 
                    theta.grid2[k], length.out = 25))
                  tile.x2 <- length.ring2 * cos(pi/2 - seq(theta.grid1[k], 
                    theta.grid2[k], length.out = 25))
                  tile.y1 <- length.ring1 * sin(pi/2 - seq(theta.grid1[k], 
                    theta.grid2[k], length.out = 25))
                  tile.y2 <- length.ring2 * sin(pi/2 - seq(theta.grid1[k], 
                    theta.grid2[k], length.out = 25))
                  polygon(c(tile.x1, rev(tile.x2)), c(tile.y1, 
                    rev(tile.y2)), col = grid.colors[k, i], border = NA)
                }
                if (grid.horizontal) {
                  segments(x0 = length.ring1 * cos(pi/2 - theta.grid1), 
                    y0 = length.ring1 * sin(pi/2 - theta.grid1), 
                    x1 = length.ring2 * cos(pi/2 - theta.grid1), 
                    y1 = length.ring2 * sin(pi/2 - theta.grid1), 
                    col = grid.col, lty = grid.lty)
                  segments(x0 = length.ring1 * cos(pi/2 - theta.grid2[length(theta.grid2)]), 
                    y0 = length.ring1 * sin(pi/2 - theta.grid2[length(theta.grid2)]), 
                    x1 = length.ring2 * cos(pi/2 - theta.grid2[length(theta.grid2)]), 
                    y1 = length.ring2 * sin(pi/2 - theta.grid2[length(theta.grid2)]), 
                    col = grid.col, lty = grid.lty)
                }
                if (grid.vertical) {
                  if (i > 1) {
                    lines((length.ring1 - length.intergr/2 + 
                      0.3 * length.intergr) * cos.tsoft, (length.ring1 - 
                      length.intergr/2 + 0.3 * length.intergr) * 
                      sin.tsoft, col = grid.col, lty = grid.lty)
                  }
                }
            }
            if (plot.type == "barplot" | plot.type == "dotplot") {
                options(warn = -1)
                if (!is.null(error.bar.inf)) {
                  arrows(x0 = length.values * cos.t, x1 = length.arrow.inf * 
                    cos.t, y0 = length.values * sin.t, y1 = length.arrow.inf * 
                    sin.t, lwd = 1, col = error.bar.col, angle = 90, 
                    length = 0.04)
                }
                if (!is.null(error.bar.sup)) {
                  arrows(x0 = length.values * cos.t, x1 = length.arrow.sup * 
                    cos.t, y0 = length.values * sin.t, y1 = length.arrow.sup * 
                    sin.t, lwd = 1, col = error.bar.col, angle = 90, 
                    length = 0.04)
                }
                options(warn = 1)
            }
            if (show.box) {
                if (tree.open.crown) {
                  lines(c(xx1, rev(xx2)), c(yy1, rev(yy2)), col = 1)
                }
                else {
                  lines(xx1, yy1, col = 1)
                  lines(xx2, yy2, col = 1)
                }
            }
            if (show.trait) {
                if (tree.open.crown) {
                  theta.trait <- theta.soft[length(theta.soft)]
                }
                else {
                  if (show.data.axis & (plot.type == "barplot" | 
                    plot.type == "dotplot")) {
                    theta.trait <- theta[length(theta)] - (360 - 
                      real.open) * (1/3) * pi/180
                  }
                  else {
                    theta.trait <- theta[length(theta)] - (360 - 
                      real.open) * (1/2) * pi/180
                  }
                }
                cos.ttrait <- cos(pi/2 - theta.trait)
                sin.ttrait <- sin(pi/2 - theta.trait)
                segments(x0 = length.ring1 * cos.ttrait, x1 = length.ring2 * 
                  cos.ttrait, y0 = length.ring1 * sin.ttrait, 
                  y1 = length.ring2 * sin.ttrait, lwd = 20, col = trait.bg.col[i])
                text(x = (length.ring1 + length.ring2)/2 * cos.ttrait, 
                  y = (length.ring1 + length.ring2)/2 * sin.ttrait, 
                  labels = trait.labels[i], col = trait.col[i], 
                  cex = trait.cex[i], font = trait.font[i], srt = ifelse(theta.trait > 
                    0 | theta.trait < -pi, (pi/2 - theta.trait) * 
                    180/pi, (-pi/2 - theta.trait) * 180/pi))
            }
        }
        if (show.tip) {
            length.tipsline <- (length.phylo + length.intergr * 
                (n.traits + 1) + length.gr * (n.traits))
            tip.xlim <- c(-1, 1)
            if (tip.adj < 0.5) 
                tip.xlim[1] <- -tip.adj/0.5
            if (tip.adj > 0.5) 
                tip.xlim[2] <- -2 * tip.adj + 2
            for (i in 1:n.tips) {
                text(x = length.tipsline * cos.t[i], y = length.tipsline * 
                  sin.t[i], labels = tip.labels[i], adj = ifelse(theta[i] > 
                  0 | theta[i] < -pi, 0, 1), col = tip.col[i], 
                  cex = tip.cex[i], font = tip.font[i], srt = ifelse(theta[i] > 
                    0 | theta[i] < -pi, (pi/2 - theta[i]) * 180/pi, 
                    (-pi/2 - theta[i]) * 180/pi))
            }
        }
    }
    par(mar = par.mar0, xpd = par.xpd0, lend = par.lend0)
    invisible()
}
