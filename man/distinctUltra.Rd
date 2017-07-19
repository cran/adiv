\name{distinctUltra}
\alias{distinctUltra}
\title{
Ultrametric Tree-based Species' Originality 
}
\description{
The function calculates two indices of species' originality that rely on the structure and branch lengths of ultrametric (phylogenetic) trees.
}
\usage{
distinctUltra(phyl, method = 1:2)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}, with ultrametric properties.
}
  \item{method}{
a numeric or a vector of numerics. Possible values are \code{1}, and \code{2}. \code{1} is for the Pavoine et al. QE-based (also named Qb) index (Pavoine et al. 2005); \code{2} is for Pavoine et al. 2H-based index (which could also be named more shortly as 2Hb) (vector that maximizes index \eqn{^2H}{2H}, Pavoine and Izsak 2014).
}
}
\value{
A data frame with species as rows and originality indices as columns.
}
\references{
Pavoine, S., Ollier, S., Dufour, A.B. (2005) Is the originality of a species measurable? \emph{Ecology Letters}, \bold{8}, 579--586.

Pavoine, S. and Izsak, J. (2014) New biodiversity measure that includes consistent interspecific and intraspecific components. \emph{Methods in Ecology and Evolution}, \bold{5}, 165--172.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}

with contributions of Stephane Dray
}
\seealso{
See Also \code{\link{distinctDis}}, \code{\link{distinctTopo}}, \code{\link{distinctTree}}
}
\examples{
data(carni70, package = "adephylo")
tre <- read.tree(text=carni70$tre)
U <- distinctUltra(tre)
U.4d <- phylo4d(tre, as.matrix(U))
parmar <- par()$mar
par(mar=rep(.1,4))
table.phylo4d(U.4d, show.node=FALSE, symbol="squares", center=FALSE, scale=FALSE, cex.symbol=2)
par(mar=parmar)
}
\keyword{models}
