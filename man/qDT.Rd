\name{qDT}
\alias{qDT}
\title{
Mean Diversity along a Phylogenetic Tree
}
\description{
Function \code{qDT} calculates the \eqn{^qD(T)} index developed by Chao et al. (2010) as the mean diversity of order \emph{q} over \emph{T} years in a phylogenetic tree. In function \code{qDT}, the index is computed over the whole tree from root to tips. It uses the formula of the \eqn{^qD(T)} index extended to non-ultrametric trees (where the distance from tip to root varies) (Chao et al. 2010, equation 4.5) . 
}
\usage{
qDT(phyl, comm, q = 2, tol = 1e-08)
}
\arguments{
  \item{phyl}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.}
  \item{comm}{a data frame or a matrix typically with communities as rows, species as columns and abundance as entry. Species should be labeled as in the phylogenetic tree where they are the tips.}
  \item{q}{a vector with nonnegative value(s) for parameter \code{q}. See details.}
  \item{tol}{numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.}
}
\value{
If only one value of \code{q} is given, a vector with the phylogenetic diversity of each community is returned.
If more than one value of \code{q} is given, a list of two objects is returned: 
\item{q}{the vector of values for \code{q};}
\item{div}{a data frame with the phylogenetic diversity of each community calculated for all values of \code{q}.}

The results of function \code{plot.qDT} are of class \code{"evodivparam"}. As such, function \code{\link{plot.evodivparam}} allows to plot the results of function \code{qDT}.
}
\references{
Chao, A., Chiu, C.-H., Jost, L. (2010) Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society London Series B}, \bold{365}, 3599--3609.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evodivparam}}
}
\examples{
\dontrun{
if(require(ape)){

data(batcomm)
phy <- read.tree(text=batcomm$tre)
ab <- batcomm$ab[, phy$tip.label]
plot(qDT(phy, ab))
plot(qDT(phy, ab, q=seq(0, 10, length=20)))

}
}
}
\keyword{models}
