\name{DP}
\alias{DP}
\title{Plot-to-plot functional or phylogenetic dissimilarity}
\description{
The function \code{DP} calculates Ricotta et al. (2020) plot-to-plot functional or phylogenetic beta uniqueness (index named \eqn{D_F} for functional data and \eqn{D_P} for phylogenetic data, calculated with equation 2 in Ricotta et al. 2020).
}
\usage{
DP(mtree, comm, height = NULL, diag = FALSE, upper = FALSE, tol = 0.001)
}
\arguments{
\item{mtree}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}. The tree must be ultrametric: equal distance from any tip to the root.}
\item{comm}{a matrix containing the relative or absolute abundance of all species in plots. Columns are species and plots are rows. Column labels (species names) should be assigned as in \code{mtree}.}
\item{height}{either \code{NULL} or a numeric. See details.}
\item{diag}{a logical value indicating whether the diagonal of the distance matrix should be printed by function \code{\link{print.dist}}.}
\item{upper}{a logical value indicating whether the upper triangle of the distance matrix should be printed by function \code{\link{print.dist}}.}
\item{tol}{a tolerance threshold. A value between -\code{tol} and \code{tol} is considered as zero. See details.}
}
\details{
Object \code{mtree} defines a tree with species as tips. If argument \code{height} is \code{NULL}, then the root of the tree will be placed at the most recent common ancestor of all species occurring in the set of plots (given in object \code{comm}). An alternative position for the root can be given by specifying the height of the tree (argument \code{height}). In that case, \code{height} must be higher than the distance between tips and the most recent common ancestor of all species.

The tolerance threshold \code{tol} is particularly important if your tree is not exactly ultrametric due to approximation problems. In that case, the distance from tip to root varies according to the tip considered, although it should not (variations are due to approximation problems). A difference smaller than \code{tol} in the  distance to root for two species will thus be considered as null.
}
\value{
The function returns a (semi-)matrix of class \code{dist} with the values of the proposed dissimilarities for each pair of plots. 
}
\references{
Ricotta, C., Laroche, F., Szeidl, L., Pavoine, S. (2020) From alpha to beta functional and phylogenetic redundancy. \emph{Methods in Ecology and Evolution}. In press.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{betaTreeUniqueness}} for beta uniqueness, \code{\link{treeUniqueness}} for alpha uniqueness
}
\examples{
\dontrun{
if(require(ape) && require(ade4)){
data(RutorGlacier)
phy <- read.tree(text=RutorGlacier$TreeNW)
plot(phy)

ab <- RutorGlacier$Abund[, phy$tip.label]

# Phylogenetic dissimilarities between plots
# (Ricotta et al. 2020)
Dp <- DP(phy, ab, tol=0.00001)
# Principal Coordinate Analysis of the dissimilarities
pcoDp <- dudi.pco(sqrt(Dp), full=TRUE)
s.class(pcoDp$li, as.factor(RutorGlacier$Fac))
}
}
}
