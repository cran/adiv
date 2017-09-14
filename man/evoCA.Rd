\name{evoCA}
\alias{evoCA}
\title{
Evolutionary Correspondence Analysis
}
\description{
The function performs evoCA (Pavoine 2016): an adaptation of the correspondence analysis (CA) to analyse the distributions of lineages among sites and, simultaneously, to analyse the phylogenetic composition of sites.
}
\usage{
evoCA(phyl, comm, scannf = TRUE, nf = 2, abundance = TRUE)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}. To ease the interpretation of factorial maps, it is advised to add labels to the nodes of the phylogeny. For example, If your phylogeny is named 'tree' and is of class \code{phylo}, you can name the nodes with the following command:
\code{tree$node.label <- paste("n", 1:tree$Nnode, sep="")}. If it is of class \code{phylo4}, you can use: \code{nodeLabels(tree) <- paste("n", 1:nNodes(tree), sep="")}. Use \code{plot(tree)} to see the result.
}
  \item{comm}{
a data frame or a matrix typically with communities (or sites, plots, etc.) as rows, species as columns and presence/absence or an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips.
}
  \item{scannf}{
a logical value indicating whether the screeplot (eigenvalues) should be displayed for choosing the number of axes to be kept.
}
  \item{nf}{
if \code{scannf} is \code{FALSE}, an integer indicating the number of axes to be kept.
}
  \item{abundance}{
a logical value, if \code{FALSE}, only the presences/absences of the \bold{evolutionary units} in sites are considered.}
}
\value{
an object of class \code{dudi} (see package ade4, \code{?dudi}). Graphical tools are associated with this class (see example section below and \code{?scatter.dudi}).
}
\references{
Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
See Also \code{\link{evoCA}}, \code{\link{evopcachord}}, \code{\link{evopcahellinger}}, \code{\link{evodiss}}
}
\examples{
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
}
\keyword{models}
