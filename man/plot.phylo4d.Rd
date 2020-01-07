\name{plot.phylo4d}
\alias{plot.phylo4d}
\alias{barp4d}
\alias{dotp4d}
\alias{gridp4d}
\title{
Plots for phylogenetic tree and data
}
\description{
Functions \code{barp4d}, \code{dotp4d} and \code{gridp4d} provide plots for phylo4d objects (i.e. phylogenetic tree and data). Function \code{plot.phylo4d} provides a general interface for all other functions. Function \code{barp4d} uses barplots of trait values along the phylogenetic tree. Function \code{dotp4d} uses dotplots of trait values along the phylogenetic tree, and function \code{gridp4d} gridplots. 
}
\usage{
\method{plot}{phylo4d}(x, trait = names(tdata(p4d)), center = TRUE, 
    scale = TRUE, plot.type = "barplot", tree.ladderize = FALSE, 
    tree.type = "phylogram", tree.ratio = NULL, 
    tree.xlim = NULL, tree.open.angle = 0, 
    tree.open.crown = TRUE, show.tip = TRUE, tip.labels = NULL, 
    tip.col = "black", tip.cex = 1, tip.font = 3, tip.adj = 0, 
    data.xlim = NULL, bar.lwd = 10, bar.col = "grey35", 
    show.data.axis = TRUE, dot.col = "black", dot.pch = 20, 
    dot.cex = 2, cell.col = topo.colors(100), 
    show.color.scale = TRUE, show.trait = TRUE, 
    trait.labels = NULL, trait.col = "black", 
    trait.cex = 1, trait.font = 1, trait.bg.col = "grey90", 
    error.bar.sup = NULL, error.bar.inf = NULL, 
    error.bar.col = 1, show.box = FALSE, grid.vertical = TRUE,
    grid.horizontal = FALSE, grid.col = "grey25", 
    grid.lty = "dashed", ...) 

barp4d(height, trait = names(tdata(height)), 
    center = TRUE, scale = TRUE, tree.ladderize = FALSE, 
    tree.type = "phylogram", tree.ratio = NULL, 
    tree.xlim = NULL, tree.open.angle = 0, 
    tree.open.crown = TRUE, show.tip = TRUE, tip.labels = NULL, 
    tip.col = "black", tip.cex = 1, tip.font = 3, tip.adj = 0, 
    data.xlim = NULL, bar.lwd = 10, bar.col = "grey35", 
    show.data.axis = TRUE, show.trait = TRUE, 
    trait.labels = NULL, trait.col = "black", trait.cex = 1,
    trait.font = 1, trait.bg.col = "grey90", 
    error.bar.sup = NULL, error.bar.inf = NULL, 
    error.bar.col = 1, show.box = FALSE, grid.vertical = TRUE, 
    grid.horizontal = FALSE, grid.col = "grey25", 
    grid.lty = "dashed", ...) 

dotp4d(p4d, trait = names(tdata(p4d)), center = TRUE, 
    scale = TRUE, tree.ladderize = FALSE, 
    tree.type = "phylogram", tree.ratio = NULL, 
    tree.xlim = NULL, tree.open.angle = 0, 
    tree.open.crown = TRUE, show.tip = TRUE, tip.labels = NULL,
    tip.col = "black", tip.cex = 1, tip.font = 3, 
    tip.adj = 0, data.xlim = NULL, show.data.axis = TRUE, 
    dot.col = "black", dot.pch = 20, dot.cex = 2, 
    show.trait = TRUE, trait.labels = NULL, trait.col = "black", 
    trait.cex = 1, trait.font = 1, trait.bg.col = "grey90", 
    error.bar.sup = NULL, error.bar.inf = NULL, 
    error.bar.col = 1, show.box = FALSE, grid.vertical = FALSE,
    grid.horizontal = TRUE, grid.col = "grey25", 
    grid.lty = "dashed", ...) 

gridp4d(p4d, trait = names(tdata(p4d)), center = TRUE, 
    scale = TRUE, tree.ladderize = FALSE, 
    tree.type = "phylogram", tree.ratio = NULL, 
    tree.xlim = NULL, tree.open.angle = 0, 
    tree.open.crown = TRUE, show.tip = TRUE, tip.labels = NULL, 
    tip.col = "black", tip.cex = 1, tip.font = 3, tip.adj = 0,
    cell.col = topo.colors(100), show.color.scale = TRUE, 
    show.trait = TRUE, trait.labels = NULL, trait.col = "black", 
    trait.cex = 0.7, trait.font = 1, trait.bg.col = "grey90", 
    show.box = FALSE, grid.vertical = FALSE, 
    grid.horizontal = FALSE, grid.col = "grey25", 
    grid.lty = "dashed", ...) 

}
\arguments{
\item{x, p4d, height}{a \code{phylo4d} object.}
\item{trait}{the traits in the \code{phylo4d} object to include in the plot. Can be a character vector giving the name of the traits or numbers giving the column index in the table of the data slot of the \code{p4d} or \code{height} object. Can be used to reorder the traits in the plot.}
\item{center}{a logical indicating whether traits values should be centered.}
\item{scale}{a logical indicating whether traits values should be scaled.}
\item{plot.type}{a character string specifying the type of plot for traits data. Can be \code{"barplot"}, \code{"dotplot"} or \code{"gridplot"}.}
\item{tree.ladderize}{a logical indicating whether the tree should be (right) ladderized.}
\item{tree.type}{a character string specifying the type of phylogeny to be drawn. Can be \code{"phylogram"}, \code{"cladogram"} or \code{"fan"}.}
\item{tree.ratio}{a numeric value in [0, 1] giving the proportion of width of the figure for the tree.}
\item{tree.xlim}{a numeric vector of length 2 giving the limits of the x-axis for the tree. If NULL, it is determined automatically.}
\item{tree.open.angle	}{a numeric value giving the angle in degrees left blank if \code{tree.type = "fan"}.}
\item{tree.open.crown}{a logical indicating whether the crowns should be drawn following the value of \code{tree.open.angle} (default \code{TRUE}).}
\item{show.tip}{logical indicating whether tips labels should be drawn.}
\item{tip.labels}{character vector to label the tips. If \code{NULL} the tips labels of the \code{phylo4d} object are used.}
\item{tip.col}{a vector of R colors to use for the tips labels. Recycled if necessary.}
\item{tip.cex}{a numeric vector to control character size of the tips labels. Recycled if necessary.}
\item{tip.font}{an integer vector specifying the type of font for the tips labels: 1 (plain text), 2 (bold), 3 (italic), or 4 (bold italic). Recycled if necessary.}
\item{tip.adj}{a vector of numeric in [0, 1] to control tips labels justification: 0 (left-justification), 0.5 (centering), or 1 (right-justification). Recycled if necessary.}
\item{data.xlim}{numeric vector of length 2 or matrix giving the x coordinates range for the barplots/dotplots.}
\item{bar.lwd}{a vector of numeric giving bar widths of the barplot(s). Recycled along the tips, reapeated for each trait.}
\item{bar.col}{a vector of R colors to use for the bars. Recycled along the tips, reapeated for each trait. The user can also provide a matrix for a finer tuning.}
\item{show.data.axis}{logical indicating whether barplots/dotplots axes should be drawn.}
\item{dot.col}{a vector of R colors to use for the points. Recycled along the tips, reapeated for each trait. The user can also provide a matrix for a finer tuning.}
\item{dot.pch}{a numerical vector of symbol to use for the points. Recycled along the tips, reapeated for each trait. The user can also provide a matrix for a finer tuning.}
\item{dot.cex}{a numerical vector. Character (or symbol) expansion for the points. Recycled along the tips, reapeated for each trait. The user can also provide a matrix for a finer tuning.}
\item{cell.col}{a vector of colors for \code{gridplot} cells. Easily generated by \code{\link{heat.colors}}, \code{\link{topo.colors}}, \code{\link{terrain.colors}} or other functions created with \code{link{colorRampPalette}}.}
\item{show.color.scale}{logical indicating whether color scale should be drawn.}
\item{show.trait}{logical indicating whether traits labels should be drawn.}
\item{trait.labels}{character vector to label the traits. If \code{NULL} the traits labels of the \code{phylo4d} object are used.}
\item{trait.col}{a vector of R colors to use for the traits labels. Recycled if necessary.}
\item{trait.cex}{a numeric vector to control character size of the trait labels. Recycled if necessary.}
\item{trait.font}{an integer vector specifying the type of font for the traits labels: 1 (plain text), 2 (bold), 3 (italic), or 4 (bold italic). Recycled if necessary.}
\item{trait.bg.col}{a vector of R colors to use for the background of the barplots. Recycled if necessary.}
\item{error.bar.sup}{a matrix giving the superior limit for error bars. Columns and rows names must match with traits and tips labels, respectively.}
\item{error.bar.inf}{	a matrix giving the inferior limit for error bars. Columns and rows names must match with traits and tips labels, respectively.}
\item{error.bar.col}{	a vector of R colors to draw error bars.}
\item{show.box}{a logical indicating whether a box should be drawn around the plots.}
\item{grid.vertical}{	a logical incating whether vertical lines of the grid should be drawn.}
\item{grid.horizontal}{a logical incating whether horizontal lines of the grid should be drawn.}
\item{grid.col}{a vector of R colors to use for the lines of the grid.}
\item{grid.lty}{the lines type of the grid. Possibly a vector.}
\item{\dots}{further arguments to be passed to \code{plot.phylo}.}
}
\author{
The four functions were written by Francois Keck in the package named phylosignal. Functions were there named as follows: multiplot.phylo4d, barplot.phylo4d, dotplot.phylo4d, and gridplot.phylo4d.
At the end of 2019, the package was orphaned and the functions were integrated in package \code{adiv}. The versions of the functions have been slightly modified compared to those developed by Francois Keck.
}
\examples{
if(require(ape) && require(phylobase)){
data(batcomm)
# Phylogenetic tree for bat species
phy <- read.tree(text = batcomm$tre)
# Abondance data plotted in front of the phylogeny
# F = rainforest, P = cacao plantation
# O = oldfields, C = cornfields
ab.4d <- phylo4d(phy, t(batcomm$ab))
barp4d(ab.4d, center = FALSE, scale = FALSE, 
    data.xlim = c(0, max(batcomm$ab)))
dotp4d(ab.4d, center = FALSE, scale = FALSE, 
    data.xlim = c(0, max(batcomm$ab)))
gridp4d(ab.4d, center = FALSE, scale = FALSE)
}
}
\keyword{models}
