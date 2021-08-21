\name{FPdivparam}
\alias{FPdivparam}
\alias{plot.FPdivparam}
\title{
Parametric Indices of Functional and Phylogenetic Diversity
}
\description{
Function \code{FPivparam} calculates functional or phylogenetic (FP-)diversity in communities using parametric indices \eqn{^{\alpha}K}{alphaK}, \eqn{^{\alpha}K*}{alphaK*}, \eqn{^{\alpha}Y}{alphaY} and \eqn{^{\alpha}I}{alphaI} discussed and developed in Pavoine and Ricotta (2021). 
The function \code{plot.FPdivparam} plots the results of function \code{FPdivparam}.
}
\usage{
FPdivparam(comm, disORtree, method = c("KY", "KstarI"), palpha = 2, 
    equivalent = FALSE, option = c("asymmetric", "symmetric"), 
    dmax = NULL, tol = 1e-8)

\method{plot}{FPdivparam}(x, legend = TRUE, 
    legendposi = "topright", axisLABEL = "FP-diversity", 
    type="b", col = if(is.numeric(x)) NULL 
    else sample(colors(distinct = TRUE), nrow(x$div)), 
    lty = if(is.numeric(x)) NULL else rep(1, nrow(x$div)), 
    pch = if(is.numeric(x)) NULL else rep(19, nrow(x$div)), 
    \dots)
}
\arguments{
  \item{comm}{a data frame or a matrix typically with communities as rows, species as columns and abundance as entry. Species should be labelled as in object disORtree.}
  \item{disORtree}{an object inheriting the class \code{dist}, giving the FP-dissimilarities between species, or inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}, where species are tips.}
  \item{method}{a string either \code{"KY"} for index \eqn{^{\alpha}K}{alphaK} (if \code{disORtree} is of class \code{dist}) or \eqn{^{\alpha}Y}{alphaY} if \code{disORtree} is of class \code{phylo}, \code{phylo4} or \code{hclust}; or \code{"KstarI"} for index \eqn{^{\alpha}K*}{alphaK*} (if \code{disORtree} is of class \code{dist}) or \eqn{^{\alpha}I}{alphaI} if \code{disORtree} is of class \code{phylo}, \code{phylo4} or \code{hclust}. If several values are given, only the first one is considered.}
  \item{palpha}{a vector with nonnegative value(s) for parameter \eqn{\alpha}.}
  \item{equivalent}{a logical. If TRUE, the diversity values are calculated in terms of equivalent number of species (see Pavoine and Ricotta, 2021).}
  \item{option}{a string either \code{"asymmetric"} or \code{"symmetric"}. The parameter \code{option} is only used if \code{disORtree} is of class \code{phylo}, \code{phylo4} or \code{hclust} If \code{option="symmetric"}, the distance between two tips on a tree is defined as half the sum of branch length on the smallest path that connects the two species; while if \code{option="asymmetric"} the distance between a tip \emph{i} and a tip \emph{j} on a tree is defined as the sum of branch lengths between tip \emph{i} and its most recent ancestor with tip \emph{j}. If the tree is ultrametric, the two options are equivalent.}
  \item{dmax}{a nonnegative numeric indicating the maximum possible dissimilarity between two species. \code{dmax} must be higher than or equal to the maximum observed dissimilarity between two species.}
  \item{tol}{a numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.}
  \item{x}{an object of class \code{FPdivparam} obtained with function \code{FPdivparam}.}
  \item{legend}{a logical. If TRUE a legend is given with the colour, the type of line (etc.) used to define the diversity curve of each community.}
  \item{legendposi}{a string that gives the position of the legend to be passed to function \code{legend} of the base of R.}
\item{axisLABEL}{a string to display on the main axis of the plot to designate what we are measuring. The default is \code{"FP-diversity"}.} 
  \item{type}{a string to be passed to the graphic argument \code{type} of functions \code{plot} and \code{lines} used to draw the diversity curve of each community.}
  \item{col}{vector of colours to be passed to the graphic argument \code{col} of functions \code{plot} and \code{lines} to define the colour of the diversity curve of each community.}
  \item{lty}{vector of type of line (plain, broken etc.) to be passed to the graphic argument \code{lty} of functions \code{plot} and \code{lines} used to draw the diversity curve of each community.}
  \item{pch}{type of point (open circle, close circle, square etc.) to be passed to the graphic argument \code{pch} of functions \code{plot} and \code{lines} used to draw the diversity level of each community.}
  \item{\dots}{other arguments can be added and passed to the functions \code{plot} and \code{lines} used to draw the graphic.}
}
\value{
If only one value of \code{palpha} is given, function \code{FPdivparam} returns a vector with the phylogenetic diversity of each community.
If more than one value of \code{palpha} is given, a list of two objects is returned: 
\item{palpha}{the vector of values for \code{palpha};}
\item{div}{a data frame with the phylogenetic diversity of each community calculated for all values of \code{palpha}.}

The function \code{plot.FPdivparam} returns a graphic.
}
\references{
The methodologies and scripts were developed by

Pavoine, S., Ricotta, C. (2021) On the relationships between rarity, uniqueness, distinctiveness, originality and functional/phylogenetic diversity. \emph{BiorXiv}. \doi{10.1101/2021.08.09.455640}

unifying and extending earlier work by:

Pavoine, S., Love, M., Bonsall, M. (2009) Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: application to rockfish (genus: Sebastes) in the Southern California Bight. \emph{Ecology Letters}, \bold{12}, 898--908.

and

Ricotta, C., Szeidl, L. (2006) Towards a unifying approach to diversity measures: Bridging the gap between the Shannon entropy and Rao's quadratic index. \emph{Theoretical Population Biology}, \bold{70}, 237--243.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{pIa}}, \code{\link{evodivparam}}
}
\examples{
\dontrun{
if(require(ape)){

data(batcomm)
phy2 <- read.tree(text=batcomm$tre2)
ab2 <- batcomm$ab2[, phy2$tip.label]
plot(FPdivparam(ab2, phy2))
plot(FPdivparam(ab2, phy2, palpha=seq(0, 10, length=20)))

}
}
}
\keyword{models}
