\name{distinctTree}
\alias{distinctTree}
\title{
Tree-based Species' Originality 
}
\description{
The function calculates indices of species' originality that rely on the structure and branch lengths of (phylogenetic) trees. Trees with polytomies are allowed.
}
\usage{
distinctTree(phyl, method = c("ED", "ES", "Delta*"), palpha = 0, standardized = FALSE)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}. 
}
  \item{method}{
a string or a vector of strings. Possible values are \code{"ED"}, 
\code{"ES"}, and \code{"Delta*"}. \code{"ED"} is for the evolutionary 
distinctiveness, also named fair-proportion, index (Redding 2003; Isaac et 
al. 2007); \code{"ES"} is for the Equal-Splits index (Redding and Mooers 
2006); \code{"Delta*"} is for the parametric \eqn{^{\alpha}\Delta} index of 
Pavoine and Ricotta (2021).}
  \item{palpha}{a numeric value or a numeric vector of values for parameter \eqn{\alpha} of \eqn{^{\alpha}\Delta} index by 
Pavoine and Ricotta (2021).} 
  \item{standardized}{
a logical. If \code{TRUE}, the vector of originalities is divided by its sum (transforming absolute originalities into relative originalities).
}
}
\value{
A data frame with species as rows and originality indices as columns.
}
\references{
Isaac, N.J., Turvey, S.T., Collen, B., Waterman, C., Baillie, J.E. (2007) Mammals on the EDGE: conservation priorities based on threat and phylogeny. \emph{PloS ONE}, \bold{2}, e296.

Pavoine, S., Ricotta, C. (2021) On the relationships between rarity, uniqueness, distinctiveness, originality and functional/phylogenetic diversity. \emph{BiorXiv}. \doi{10.1101/2021.08.09.455640}

Redding, D.W. (2003) \emph{Incorporating genetic distinctness and reserve occupancy into a conservation priorisation approach}. Master thesis: University of East Anglia, Norwich.

Redding, D.W., Mooers, A.O. (2006) Incorporating evolutionary measures into conservation prioritization. \emph{Conservation Biology}, \bold{20}, 1670--1678.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{distinctDis}}, \code{\link{distinctTopo}}, \code{\link{distinctUltra}}
}
\examples{
\dontrun{
if(require(ape) && require(adephylo) && require(phylobase)){
data(carni70, package = "adephylo")
tre <- read.tree(text=carni70$tre)
U <- distinctTree(tre, standardize = TRUE)
U.4d <- phylo4d(tre, as.matrix(U))
dotp4d(U.4d, center = FALSE, scale = FALSE)
barp4d(U.4d, center = FALSE, scale = FALSE)
gridp4d(U.4d, center = FALSE, scale = FALSE)

parmar <- par()$mar
par(mar=rep(.1,4))
table.phylo4d(U.4d, show.node=FALSE, symbol="squares", center=FALSE, scale=FALSE, cex.symbol=2)
par(mar=parmar)
}
}
}
\keyword{models}
