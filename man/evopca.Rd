\name{evopca}
\alias{evopcachord}
\alias{evopcahellinger}
\title{
Evolutionary Principal Component analyses
}
\description{
The functions \code{evopcachord} and \code{evopcahellinger} perform two adaptations of Principal Component Analysis (PCA) for the analysis of phylogenetic diversity patterns across species communities: the evolutionary PCA based on Chord distance (evoPCAChord) and the evolutionary PCA based on Hellinger distance (evoPCAHellinger) (Pavoine 2016).
}
\usage{
evopcachord(phyl, comm, option = c("centred", "decentred"), 
w = c("evoab", "even", "speciesab"), scannf = TRUE, 
nf = 2, abundance = TRUE)

evopcahellinger(phyl, comm, option = c("centred", "decentred"),
w = c("evoab", "even", "speciesab"), scannf = TRUE, 
nf = 2, abundance = TRUE)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}. To ease the interpretation of factorial maps, it is advised to use a rooted phylogenetic tree and to add labels to the nodes of the phylogeny. For example, If your phylogeny is named 'tree' and is of class \code{phylo}, you can root the phylogeny using \code{tree$root.edge <- 0} and you can name the nodes with the following command:
\code{tree$node.label <- paste("n", 1:tree$Nnode, sep="")}. Use \code{plot(tree, show.node.label=TRUE)} to see the result. If it is rooted and of class \code{phylo4}, you can use: \code{nodeLabels(tree) <- paste("n", 1:nNodes(tree), sep="")}. 
}
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and presence/absence (1/0) or an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips.
}
  \item{option}{
a string: either \code{"centered"} or \code{"decentered"}. If several strings are given, only the first one is used. See details.}
  \item{w}{
a string: either \code{"evoab"}, \code{"even"}, or \code{"speciesab"}. If several strings are given, only the first one is used. See details.
}
  \item{scannf}{
a logical value indicating whether the screeplot (eigenvalues) should be displayed for choosing the number of axes to be kept.
}
  \item{nf}{
if \code{scannf} is \code{FALSE}, an integer indicating the number of axes to be kept.
}
  \item{abundance}{
a logical value. If \code{FALSE}, only the presences/absences of the \bold{evolutionary units} in communities are considered.}
}
\details{
Let \eqn{L_k} be the length of branch \emph{k} in the phylogenetic tree (out of \emph{K} branches); \eqn{a_{jk}}{a_jk} the sum of abundances, in community \emph{j}, for all species descending from branch \emph{k}; \eqn{w_j} a positive weight attributed to community \emph{j} (the definition for \eqn{w_j} is flexible with the only requirement that \eqn{\sum_{j=1}^m w_j=1}{sum_j w_j=1}); 
\eqn{a_{j+}=\sum_{k=1}^K L_k a_{jk}}{a_j+=sum_k L_k a_jk}; \eqn{a_{+k}=\sum_{j=1}^m a_{jk}}{a_+k=sum_j a_jk}; \eqn{a_{++}=\sum_{j=1}^m\sum_{k=1}^K L_k a_{jk}}{a_++=sum_j sum_k L_k a_jk}. 

The weights of the communities (argument \code{w}) can be \code{"even"} (even weights, i.e. relative abundances
are considered for evolutionary units), \code{"evoab"} (proportional to the summed abundances of all evolutionary units, i.e. absolute abundances are considered for evolutionary units), or \code{"speciesab"} (proportional to the summed abundances of all species). Note that if the phylogenetic tree is ultrametric (the distance from any species to the root is constant), then options \code{"evoab"} and \code{"speciesab"} are equivalent.

In \code{evopcahellinger}, as recommended by Rao (1995), the vector used to centre matrix \eqn{(\sqrt{a_{jk}/a_{j+}})_{j,k}}{(\sqrt{a_jk/a_j+})_j,k} in PCA can be defined as
\deqn{\left(\sum_{j=1}^m w_j \sqrt{a_{jk}/a_{j+}}\right)_k}{(sum_j w_j \sqrt{a_jk/a_j+})_k} 
(ordinary weighted mean, option \code{"centered"}) 
or as 
\deqn{\left(\sqrt{(\sum_{j=1}^m w_j a_{jk})/(\sum_{k=1}^K L_k \sum_{j=1}^m w_j a_{jk})}\right)_k}{(\sqrt{(sum_j w_j a_jk)/(sum_k L_k sum_j w_j a_jk)})_k} 
(option \code{"decentered"}); see Pavoine (2016) for an introduction of all ordination approaches.

Similarly, in \code{evopcachord}, the vector used to centre matrix \eqn{(a_{jk}/\sqrt{\sum_{k=1}^K a_{jk}^2})_{j,k}}{(a_jk/\sqrt{sum_k a_jk^2})_j,k} can be defined as 
\deqn{\left(\sum_{j=1}^m w_j a_{jk}/\sqrt{\sum_{k=1}^K a_{jk}^2}\right)_k}{(sum_j w_j a_jk/\sqrt{sum_k a_jk^2})_k} 
(ordinary weighted mean, option \code{"centered"}) 
or as 
\deqn{\left((\sum_{j=1}^m w_j a_{jk})/\sqrt{\sum_{k=1}^K L_k (\sum_{j=1}^m w_j a_{jk})^2}\right)_k}{((sum_j w_j a_jk)/\sqrt{sum_k L_k (sum_j w_j a_jk)^2})_k} 
(option \code{"decentered"}).
}
\value{
\code{evopcachord} and \code{evopcahellinger} both return an object of class \code{evopca} and of class \code{dudi} (see package ade4). Graphical tools are associated with class \code{dudi} in packages ade4 and adegraphics (see section "Examples" below).

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

If X is an object of class \code{evopca}, then \code{attributes(X)$phy} contains the phylogenetic tree (of class \code{phylo}) with names for internal nodes. 
}
\references{
Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732.

Rao, C.R. (1995) A review of canonical coordinates and an alternative to correspondence analysis using Hellinger distance.
\emph{Questiio}, \bold{19}, 23--63. 
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evoCA}}, \code{\link{evoNSCA}}, \code{\link{evodiss}}
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

evopcaHbat <- evopcahellinger(phy, ab, scan=FALSE, nf=3)
dist(evopcaHbat$li)
evodiss(phy, ab, "Hellinger")

evopcaHbat$eig/sum(evopcaHbat$eig)
s.label(evopcaHbat$li)
s.label(evopcaHbat$co)
s.arrow(evopcaHbat$co)
inertia.dudi(evopcaHbat, row=TRUE)$row.abs
inertia.dudi(evopcaHbat, col=TRUE)$col.abs

evopcaCbat <- evopcachord(phy, ab, scan=FALSE, nf=3)
dist(evopcaCbat$li)
evodiss(phy, ab, "Chord")

evopcaCbat$eig/sum(evopcaCbat$eig)
s.label(evopcaCbat$li)
s.label(evopcaCbat$co)
s.arrow(evopcaCbat$co)
inertia.dudi(evopcaCbat, row=TRUE)$row.abs
inertia.dudi(evopcaCbat, col=TRUE)$col.abs

adegpar("plabels.optim" = O)
}
}
}
\keyword{models}
