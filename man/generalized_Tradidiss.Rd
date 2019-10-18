\name{generalized_Tradidiss}
\alias{generalized_Tradidiss}
\title{
Plot-to-plot dissimilarity taking account of functional dissimilarities between species
}
\description{
Given a matrix of \emph{S} species' relative or absolute abundance values in \emph{N} plots, together with an \emph{S} x \emph{S} (functional) dissimilarity matrix, the function \code{generalized_Tradidiss} calculates a semimatrix with the values of a plot-to-plot dissimilarity index, as proposed in Pavoine and Ricotta (2019).
}
\usage{
generalized_Tradidiss(comm, dis, method = c("GC", "MS", "PE"),
    abundance = c("relative", "absolute", "none"), 
    weights = c("uneven", "even"), tol = 1e-8)
}
\arguments{
  \item{comm}{
a data frame typically with communities as rows, species as columns and an index of abundance as entries. Species must be labeled as in the object \code{dis}.
}
  \item{dis}{
an object of class \code{dist} or a matrix containing pair-wise (functional or phylogenetic) dissimilarities between species [obtained by functions like vegdist in package vegan (Oksanen et al. 2013), gowdis in package FD (Laliberte and Shipley 2011), or dist.ktab in package ade4 for functional dissimilarities (Dray et al. 2007), or functions like cophenetic.phylo in package ape (Paradis et al. 2004) or distTips in package adephylo (Jombart and Dray 2010) for phylogenetic dissimilarities]. If the dissimilarities are outside the range 0-1, a warning message is displayed and each dissimilarity is divided by the maximum over all pairwise dissimilarities.
}
  \item{method}{
one of the following strings: \code{"GC"}, \code{"MS"}, \code{"PE"}. See Details.
}
  \item{abundance}{
a string with three possible values: "relative" for the use of relative species abundance, "absolute" for the use of absolute species abundance, and "none" for the use of presence/absence data (1/0). 
}
  \item{weights}{
a string. Two types of weights are available in the function: \code{"uneven"} (Eq. 5 in Pavoine and Ricotta (2009)) or \code{"even"} (1/S, where S is the number of species in the two compared plots).}
  \item{tol}{
numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.
}
}
\details{
The plot-to-plot dissimilarity coefficients used in this function are as follows:
\code{"GC"}: Equation 6 in Pavoine and Ricotta (2019)
\code{"MS"}: Equation 8 in Pavoine and Ricotta (2019)
\code{"PE"}: Equations 9 and 10 in Pavoine and Ricotta (2019)
}
\value{
The function returns an object of class \code{"dist"} with the values of the proposed dissimilarities for each pair of plots. 
}
\references{
Pavoine, S. and Ricotta, C. (2019) Measuring functional dissimilarity among plots: adapting old methods to new questions. \emph{Ecological Indicators}, \bold{97}, 67--72.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
\dontrun{
if(require(ade4) && require(adephylo) && require(ape)){

data(birdData)
phy <- read.tree(text=birdData$tre)
phydis <- sqrt(distTips(phy, method="nNodes")+1)
fau <- birdData$fau[1:6, phy$tip.label]
disGC <- generalized_Tradidiss(fau, phydis, method="GC")
disGC

### The second example is a bit TIME CONSUMING
data(mafragh)
namspe <- rownames(mafragh$traits[[1]])
M <- mafragh$flo
colnames(M) <- namspe
Bin <- prep.binary(mafragh$traits$tabBinary, c(3, 4))
distraits <- dist.ktab(ktab.list.df(list(mafragh$traits$tabOrdinal[,2:3], Bin)), 
    c("O","B"), scan=FALSE)
disGC <- generalized_Tradidiss(M, distraits, method="GC")
pcoGC <- dudi.pco(as.dist(cailliez(disGC)), full=TRUE)
s.value(mafragh$xy, pcoGC$li[,1])
}
}
}
\keyword{models}
