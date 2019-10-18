\name{evoNSCA}
\alias{evoNSCA}
\title{
Evolutionary Non-Symmetric Correspondence Analysis
}
\description{
The function performs evoNSCA (Pavoine 2016): an adaptation of the non-symmetric correspondence analysis (NSCA) (see e.g. Kroonenberg and Lombardo 1999) to analyse the distributions of lineages among sites and, simultaneously, to analyse the phylogenetic composition of sites.
}
\usage{
evoNSCA(phyl, comm, scannf = TRUE, nf = 2, abundance = TRUE)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}. To ease the interpretation of factorial maps, it is advised to use a rooted phylogenetic tree and to add labels to the nodes of the phylogeny. For example, If your phylogeny is named 'tree' and is of class \code{phylo}, you can root the phylogeny using \code{tree$root.edge <- 0} and you can name the nodes with the following command:
\code{tree$node.label <- paste("n", 1:tree$Nnode, sep="")}. Use \code{plot(tree, show.node.label=TRUE)} to see the result. If it is rooted and of class \code{phylo4}, you can use: \code{nodeLabels(tree) <- paste("n", 1:nNodes(tree), sep="")}.
}
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and presence/absence (1/0) or an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips.
}
  \item{scannf}{
a logical value indicating whether the screeplot (eigenvalues) should be displayed for choosing the number of axes to be kept.
}
  \item{nf}{
if scannf is \code{FALSE}, an integer indicating the number of axes to be kept.
}
  \item{abundance}{
a logical value, if \code{FALSE}, only the presences/absences of the \bold{evolutionary units} in sites are considered.}
}
\value{
\code{evoNSCA} returns an object of class \code{evoNSCA} and of class \code{dudi} (see package ade4). Graphical tools are associated with class \code{dudi} in packages ade4 and adegraphics (see section "Examples" below).

The returned object contains the following components:

\item{tab}{a data frame with n rows and p columns, with communities as rows and nodes of the phylogeny as columns. Internal data frame used by the algorithm;}
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

If X is an object of class \code{evoNSCA}, then \code{attributes(X)$phy} contains the phylogenetic tree (of class \code{phylo}) with names for internal nodes. 
}
\references{
Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 

Kroonenberg PM, Lombardo R (1999) Nonsymmetric correspondence analysis: a tool for analysing contingency tables with a dependence structure. \emph{Multivariate Behavioral Research}, \bold{34}, 367--396.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evoCA}}, \code{\link{evopcachord}}, \code{\link{evopcahellinger}}, \code{\link{evodiss}}
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
evoNSCAbat <- evoNSCA(phy, ab, scan=FALSE, nf=2)
evoNSCAbat$eig/sum(evoNSCAbat$eig)
s.label(evoNSCAbat$li)
s.label(evoNSCAbat$co)
s.arrow(evoNSCAbat$co)
inertia.dudi(evoNSCAbat, row=TRUE)$row.abs
inertia.dudi(evoNSCAbat, col=TRUE)$col.abs

evoNSCAbat <- evoNSCA(phy, ab, scan=FALSE, nf=3) ## All axes are now retained
# The Euclidean (canonical) distances among habitat points on the evoNSCA space is
dist(evoNSCAbat$li)
# which is equal to evoDprofile (see function evodiss):
evodiss(phy, ab, "LG")

adegpar("plabels.optim" = O)
}
}
}
\keyword{models}
