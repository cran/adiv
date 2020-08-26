\name{evoCA}
\alias{evoCA}
\alias{plot.evoCA}
\title{
Evolutionary Correspondence Analysis
}
\description{
The function \code{evoCA} performs the evolutionary correspondence analysis (evoCA) (Pavoine 2016): an adaptation of the correspondence analysis (CA) to analyse the distributions of lineages among sites and, simultaneously, to analyse the phylogenetic composition of sites. The function \code{plot.evoCA} displays the phylogeny on the factorial maps of evoCA.
}
\usage{
evoCA(phyl, comm, scannf = TRUE, nf = 2, abundance = TRUE)

\method{plot}{evoCA}(x, xaxis = 1, yaxis = 2, graph = FALSE, \dots)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}. To ease the interpretation of factorial maps, it is advised to use a rooted phylogenetic tree and to add labels to the nodes of the phylogeny. For example, If you phylogeny is named 'tree' and is of class \code{phylo}, you can root the phylogeny using \code{tree$root.edge <- 0} and you can name the nodes with the following command:
\code{tree$node.label <- paste("n", 1:tree$Nnode, sep="")}. Use \code{plot(tree, show.node.label=TRUE)} to see the result. If it is rooted and of class \code{phylo4}, you can use: \code{nodeLabels(tree) <- paste("n", 1:nNodes(tree), sep="")}.
}
  \item{comm}{
a data frame or a matrix typically with communities (or sites, plots, etc.) as rows, species as columns and presence/absence (1/0) or an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips.
}
  \item{scannf}{
a logical value indicating whether the screeplot (eigenvalues) should be displayed for choosing the number of axes to be kept.
}
  \item{nf}{
if \code{scannf} is \code{FALSE}, an integer indicating the number of axes to be kept.
}
  \item{abundance}{
a logical value, if \code{FALSE}, only the presences/absences of the \bold{evolutionary units} in sites are considered.
}
  \item{x}{
an object of class \code{evoCA}}
  \item{xaxis}{
the number of the evoCA axis chosen for the x-axis of the 3d plot.
}
  \item{yaxis}{
the number of the evoCA axis chosen for the y-axis of the 3d plot.
}
  \item{graph}{
a logical value. \code{plot.evoCA} depends on package \code{rgl}. This package will only work if you have graphics support on your computer. If you do have graphics support, choose graph = TRUE and you will be able to visualize the 3d plot with coordinates for species and phylogenetic nodes.
}
  \item{\dots}{
other arguments can be added and passed to the function \code{plot3d} of package rgl.
}
}
\value{
\code{evoCA} returns an object of class \code{evoCA} and of class \code{dudi} (see package ade4, \code{?dudi}). Graphical tools are associated with class \code{dudi} in packages ade4 and adegraphics (see example section below and \code{?scatter.dudi}).

The returned object contains the following components:

\item{tab}{a data frame with n rows and p columns, with communities as rows and nodes of the phylogeny as columns; the entries of the data frame evaluate the degree of dependence (values that depart from zero)/independence(close-to-zero values) between the occurrence in a community and the position in the phylogeny;}
\item{cw}{weights attributed to the nodes of the phylogeny, a vector with p components;}
\item{lw}{weights attributed to the communities, a vector with n components;}
\item{eig}{vector of eigenvalues;}
\item{rank}{integer, number of axes;}
\item{nf}{integer, number of kept axes;}
\item{c1}{normed scores for the nodes of the phylogeny, data frame with p rows and nf columns;}
\item{l1}{normed scores for the communities, data frame with n rows and nf columns;}
\item{co}{scores for the nodes of the phylogeny, data frame with p rows and nf columns;}
\item{li}{scores for the communities, data frame with n rows and nf columns;}
\item{call}{the original call.}

If X is an object of class \code{evoCA}, then \code{attributes(X)$phy} contains the phylogenetic tree (of class \code{phylo}) with names for internal nodes. 

\code{plot.evoCA} returns a dynamics 3-dimensional plot
}
\references{
Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evoNSCA}}, \code{\link{evopcachord}}, \code{\link{evopcahellinger}}, \code{\link{evodiss}}
}
\examples{
\dontrun{
if(require(ade4) && require(ape) && require(adegraphics)){
O <- adegpar()$plabels$optim
adegpar("plabels.optim" = TRUE)

data(batcomm)
ab <- batcomm$ab
phy <- read.tree(text=batcomm$tre)
plot(phy, show.node=TRUE)
evoCAbat <- evoCA(phy, ab, scan=FALSE, nf=2)
evoCAbat$eig/sum(evoCAbat$eig)
s.label(evoCAbat$li)
s.label(evoCAbat$co)
s.arrow(evoCAbat$co)
inertia.dudi(evoCAbat, row=TRUE)$row.abs
inertia.dudi(evoCAbat, col=TRUE)$col.abs

evoCAbat <- evoCA(phy, ab, scan=FALSE, nf=3) ## All axes are now retained
# The Euclidean (canonical) distances among habitat points on the evoCA space is
dist(evoCAbat$li)
# which is equal to evoDchi2:
evodiss(phy, ab, "chi2")

# Run the following commands only 
# if you have graphics support to use 
# the rgl package on your computer:
plot(evoCAbat, graph = TRUE)
# With argument "graph = TRUE"
# you will be able to vizualize the 3d plot.

# See also,
plot(evoCAbat, xaxis=1, yaxis=3, graph = TRUE)
adegpar("plabels.optim" = O)

}
}
}
\keyword{models}
