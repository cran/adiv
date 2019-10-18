\name{treeUniqueness}
\alias{treeUniqueness}
\title{
Community-level phylogenetic (or functional) redundancy
}
\description{
The function \code{treeUniqueness} calculates community-level phylogenetic (or tree-based) redundancy taking into account the branching pattern of the underlying phylogenetic tree (or any other tree, like a functional dendrogram).
}
\usage{
treeUniqueness(phyl, comm, index = 
    c("richness", "GiniSimpson", "Shannon"), tol = 0.001)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.
}
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and presence/absence (1/0) or an index of abundance as entries. Species should be labeled as in the tree named \code{phyl} where they are the tips.
}
  \item{index}{
a string. \code{treeUniqueness} computes tree-based (phylogenetic) species richness ("richness") or the tree-based analogue of the Shannon (shannon) or Gini-Simpson diversity ("GiniSimpson"); "richness" is assumed as the default value.
}
  \item{tol}{
a numeric. A value between -tol and tol is considered as zero. See details.
}
}
\details{
The tolerance threshold \code{tol} is particularly important if your tree is not exactly ultrametric due to approximation problems. In that case, the distance from tip to root varies according to the tip considered, although it should not (variations are due to approximation problems). A difference smaller than tol in the distance to root for two species will thus be considered as null.
}
\value{
An object of class \code{data.frame} is returned containing the following statistics:
\item{DK}{the present-day diversity of all plots in the data frame.}
\item{DP}{the tree-based (phylogenetic) diversity \eqn{D_P} of all plots in the data frame obtained by averaging the diversities \eqn{D_k} over the corresponding tree periods after rescaling the tip-to-root length of the tree to unit (for details, see Ricotta et al. 2018).}
\item{U}{tree-based (phylogenetic) uniqueness calculated as the ratio between tree-based (phylogenetic) diversity and present-day diversity \eqn{D_P/D_K}.}
\item{R}{tree-based (phylogenetic) redundancy, calculated as \eqn{1-U}.}
}
\references{
Ricotta, C., Bacaro, G., Caccianiga, M., Cerabolini, B.E.L., Pavoine, S. (2018) A new method for quantifying the phylogenetic redundancy of biological communities. \emph{Oecologia}, \bold{186}, 339--346.
}
\author{
Ricotta et al. (2018) with modifications by Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{uniqueness}}
}
\examples{
\dontrun{
if(require(ape)){
data(rockfish)
phy <- read.tree(text = rockfish$tre)
R <- treeUniqueness(phy, rockfish$fau, index = "Shannon")
}
}
}
\keyword{models}
