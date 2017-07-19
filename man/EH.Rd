\name{EH}
\alias{EH}
\title{
Faith Phylogenetic Diversity
}
\description{
Computes the sum of branch lengths on a phylogenetic tree.
}
\usage{
EH(phyl, select = NULL)
}
\arguments{
  \item{phyl}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}.}
  \item{select}{a vector containing the numbers of the leaves (species) which must be considered in the computation of Phylogenetic Diversity (PD) (or merely sum of branch lengths on the tree). This parameter allows the calculation of PD for a subset of species (including the branch between the subtree and the most ancient node of the full tree).}
}
\value{
Function \code{EH} returns a real value
}
\references{
Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity. \emph{Biological Conservation}, \bold{61}, 1--10.

Nee, S. and May, R.M. (1997) Extinction and the loss of evolutionary history. \emph{Science}, \bold{278}, 692--694. 
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
with contributions of Stephane Dray
}
\seealso{
See Also \code{\link{optimEH}}, \code{\link{randEH}}
}
\examples{
data(carni70, package = "adephylo")
tre <- read.tree(text = carni70$tre)
EH(tre)
EH(tre, select=c("Mustela.nigripes", "Mustela.frenata", "Puma.concolor"))
EH(tre, select=c(1,68,70))
}
\keyword{models}
