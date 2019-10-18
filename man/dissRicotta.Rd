\name{dissRicotta}
\alias{dissRicotta}
\title{
Plot-to-Plot Functional or Phylogenetic Dissimilarity 
}
\description{
The function calculates plot-to-plot functional or phylogenetic dissimilarity based on index \eqn{D_{AB}}{D_AB} in Ricotta et al. (2015).
}
\usage{
dissRicotta(comm, dis)
}
\arguments{
  \item{comm}{
a matrix of the relative or absolute abundance of species in communities. Columns are species and communities are rows. Column labels (species names) should be assigned as in \code{dis}.}
  \item{dis}{
a matrix of (functional or phylogenetic) dissimilarities rescaled in the range [0, 1] or an object of class \code{dist} [obtained by functions like \code{\link{dsimFun}} in this package adiv, \code{vegdist} in package vegan, \code{gowdis} in package FD, or \code{dist.ktab} in package ade4 for functional dissimilarities, or functions like \code{\link{dsimTree}} in this package adiv, \code{cophenetic.phylo} in package ape or \code{distTips} in package adephylo for phylogenetic dissimilarities; See also function \code{\link{dsimTaxo}} for taxonomic data].}
}
\value{
The function returns a semi-matrix of class \code{dist} with the values of the proposed dissimilarities for each pair of plots. Note that dissimilarities among species need first to be rescaled in the range [0, 1]. If the dissimilarities are outside the range 0-1 (as it is usually the case in phylogenetic studies for instance), a warning message is displayed and all dissimilarities are divided by the maximum observed dissimilarity.}
\references{
Ricotta, C., Bacaro, G., Pavoine, S. (2015) A cautionary note on some phylogenetic dissimilarity measures. \emph{Journal of Plant Ecology}, \bold{8}, 12--16.}
\author{
Giovanni Bacaro
and
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
\dontrun{
if(require(ape)){
# Phylogenetic tree
s<-"test(((v:20,w:20):10,(x:20,y:20):10):15,z:45):5;"
plot(test <- read.tree(text=s))

# Phylogenetic distances among species
tdist <- cophenetic(test)/100

# Matrix of abundances of the species in four communities; 
# communities A and C are identical;
# communities B and D are identical;
comm <- t(data.frame(A = rep(0.2, 5), 
B = c(0.1, 0.2, 0.2, 0, 0.5), C = rep(0.2, 5), 
D = c(0.1, 0.2, 0.2, 0, 0.5), row.names = letters[22:26]))

# Index DAB
dissRicotta(comm, tdist)
}
}
}
\keyword{models}
