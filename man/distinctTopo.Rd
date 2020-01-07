\name{distinctTopo}
\alias{distinctTopo}
\title{
Topology-based Species' Originality 
}
\description{
The function calculates three indices of species' originality that rely on the topology of (phylogenetic) trees. Trees with polytomies are allowed.
}
\usage{
distinctTopo(phyl, method = c("VW","M","A","full"), standardized = FALSE)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.
}
  \item{method}{
a character or string or a vector of characters/strings. Possible values are \code{"VW"}, \code{"M"}, \code{"A"} and \code{"full"}. \code{"VW"} is for Vane-Wright et al. index (Vane-Wright et al. 1991); \code{"M"} is for May index (May 1990); \code{"A"} is for Pavoine et al. index of originality derived from Abouheif (1999) (Pavoine et al. 2008); and "full" returns all indices.
}
  \item{standardized}{
a logical. If \code{FALSE}, the minimum score is scaled to 1 as in May (1990). If \code{TRUE}, the vector of originalities is divided by its sum (transforming absolute originalities into relative originalities).
}
}
\value{
A data frame with species as rows and originality indices as columns.
}
\references{
May, R.M. (1990) Taxonomy as destiny. \emph{Nature}, \bold{347}, 129--130.

Vane-Wright, R.I., Humphries, C.J., Williams, P.H. (1991) What to protect? Systematics and the agony of choice. \emph{Biological Conservation}, \bold{55}, 235--254.

Abouheif, E. (1999) A method for testing the assumption of phylogenetic independence in comparative data. \emph{Evolutionary Ecology Research}, \bold{1}, 895--909.

Pavoine, S., Ollier, S., Pontier, D., Chessel, D. (2008) Testing for phylogenetic signal in phenotypic traits: new matrices of phylogenetic proximities. \emph{Theoretical Population Biology}, \bold{73}, 79--91.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
with contributions of Stephane Dray
}
\seealso{
\code{\link{distinctDis}}, \code{\link{distinctTree}}, \code{\link{distinctUltra}}
}
\examples{
\dontrun{
if(require(ape) && require(adephylo) && require(phylobase)){
data(carni70, package = "adephylo")
tre <- read.tree(text=carni70$tre)
U <- distinctTopo(tre, standardized = TRUE)
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
