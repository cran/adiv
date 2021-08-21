\name{distinctAb}
\alias{distinctAb}
\title{
Abundance-based measures of species' rarity, functional or phylogenetic distinctiveness and functional or phylogenetic effective originality 
}
\description{
The function calculates parametric indices of species' rarity and functional or phylogenetic distinctiveness and effective originality according to Pavoine and Ricotta (2021).
}
\usage{
distinctAb(comm, disORtree, method = c("Q", "KY", "KstarI"), palpha = 2, 
    option = c("asymmetric", "symmetric"), tol = 1e-10)
}
\arguments{
  \item{comm}{a data frame or a matrix typically with communities as rows, species as columns and abundance as entry. Species should be labelled as in object \code{disORtree}.}
  \item{disORtree}{an object inheriting the class \code{dist}, giving the FP-dissimilarities between species, or inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}, where species are tips.}
  \item{method}{a string either \code{"Q"} for the quadratic entropy, \code{"KY"} for index \eqn{^{\alpha}K}{alphaK} (if \code{disORtree} is of class \code{dist}) or \eqn{^{\alpha}Y}{alphaY} if \code{disORtree} is of class \code{phylo}, \code{phylo4} or \code{hclust}; or \code{"KstarI"} for index \eqn{^{\alpha}K*}{alphaK*} (if \code{disORtree} is of class \code{dist}) or \eqn{^{\alpha}I}{alphaI} if \code{disORtree} is of class \code{phylo}, \code{phylo4} or \code{hclust}. If several values are given, only the first one is considered.}
  \item{palpha}{a numeric, nonnegative value for parameter \eqn{\alpha}{alpha}. \code{palpha} is ignored if \code{method = "Q"}.}
  \item{option}{a string either \code{"asymmetric"}, or\code{"symmetric"} The parameter \code{option} is only used if \code{disORtree} is of class \code{phylo}, \code{phylo4} or \code{hclust} If \code{option="symmetric"}, the distance between two tips on a tree is defined as half the sum of branch length on the smallest path that connects the two species; while if \code{option="asymmetric"} the distance between a tip \emph{i} and a tip \emph{j} on a tree is defined as the sum of branch lengths between tip \emph{i} and its most recent ancestor with tip \emph{j}. If the tree is ultrametric, the two options are equivalent.}
  \item{tol}{numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.}
}
\value{
If \code{palpha} <= 1, then, the function returns a list of four objects of class \code{data.frame} (with communities as rows and species as columns):
\code{TotContr} provides for each species in each community the effective originality multiplied by the relative abundance;
\code{EffOriPres} provides for each species present in each community the effective originality; NA for absent species;
\code{DistinctPres} provides for each species present in each community the distinctiveness; NA for absent species;
\code{Rarity} provides for each species in each community the rarity (maximum rarity for absent species).

Else, the function returns a list with the three objects of class \code{data.frame} presented above and two additional objects also of class \code{data.frame}:
\code{EffOriAll} provides for each species its effective originality compared with the composition of each community; even absent species have a value considering they have zero abundance so maximum rarity and considering their functional dissimilarity or phylogenetic dissimilarity with all species present in each community;
\code{DistinctAll} provides for each species its effective originality compared with the composition of each community; even absent species have a value considering they have zero abundance so maximum rarity and considering their functional dissimilarity or phylogenetic dissimilarity with all species present in each community.
}
\references{
Pavoine, S., Ricotta, C. (2021) On the relationships between rarity, uniqueness, distinctiveness, originality and functional/phylogenetic diversity. \emph{BiorXiv}. \doi{10.1101/2021.08.09.455640}
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{distinctDis}}, \code{\link{distinctTopo}}, \code{\link{distinctTree}}, \code{\link{distinctUltra}}
}
\examples{
\dontrun{
if(require(ape) && require(adephylo) && require(phylobase)){
data(batcomm)
phy <- read.tree(text=batcomm$tre2)
disAb <- distinctAb(batcomm$ab2, phy, method="Q") 
U.4d <- phylo4d(phy, t(disAb[[2]]))
dotp4d(U.4d, center = FALSE, scale = FALSE, data.xlim = range(t(disAb[[2]]), na.rm=TRUE), 
    tree.ratio = 0.20, dot.cex=2)
title("Effective originality associated with quadratic entropy (diversity)")
}
}
}
\keyword{models}

