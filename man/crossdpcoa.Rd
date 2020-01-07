\name{crossdpcoa}
\alias{crossdpcoa_maineffect}
\alias{crossdpcoa_version1}
\alias{crossdpcoa_version2}
\title{
Crossed-DPCoA
}
\description{
Crossed-DPCoA typically analyzes the phylogenetic or functional compositions of communities according to two factors affecting the communities (e.g. space and time; habitat and region)

The function \code{crossdpcoa_maineffect} obtains the space of DPCoA (the double principal coordinate analysis) where species are placed according to their functional traits or phylogeny and communities are placed at the center of their species. Next, levels of each factor are placed at the center of their communities. The function \code{crossdpcoa_maineffect} determines the principal axes of the positions of the levels of one of the factors in this space and projects species' points on
these principal axes. The main effect of the factor named \code{facA} is analysed by this process (Pavoine et al. 2013).

The function \code{crossdpcoa_version1} performs version 1 of the crossed DPCoA in Pavoine et al. (2013) where the effect of factor \code{facA} on the diversity of communities, given factor \code{facB}, is analysed.

The function \code{crossdpcoa_version2} performs version 2 of the crossed DPCoA in Pavoine et al. (2013) where the effect of factor \code{facA} on the diversity of communities, given factor \code{facB}, is also analysed.
}
\usage{
crossdpcoa_maineffect(df, facA, facB, dis = NULL, 
scannf = TRUE, nf = 2, w = c("classic", "independence"), 
tol = 1e-07)

crossdpcoa_version1(df, facA, facB, dis = NULL, 
scannf = TRUE, nf = 2, w = c("classic", "independence"), 
tol = 1e-07)

crossdpcoa_version2(df, facA, facB, dis = NULL, 
scannf = TRUE, nf = 2, w = c("classic", "independence"), 
tol = 1e-07)
}
\arguments{
  \item{df}{
a data frame or a matrix of 0/1 or nonnegative values. As an exemple, I consider below a communities x species data frame or matrix with abundances as entries.
}
  \item{facA}{
a factor with the same length as the number of rows (communities) in df.}
  \item{facB}{
another factor with the same length as the number of rows (communities) in df.}
  \item{dis}{
an object of class \code{dist} that contains the distances among species (e.g. functional or phylogenetic distances). If \code{NULL} equidistances are used among species. The distances must have Euclidean properties. Distances are integrated in the Euclidean Diversity Index (Champely and Chessel 2002), which corresponds to a particular formulation of Rao (1982) quadratic entropy. For example the diversity within a community \emph{i} is \eqn{EDI=\sum_{k=1}^S\sum_{l=1}^S p_{k|i}p_{k|j}\frac{d_{kl}^2}{2}}{EDI=sum_k,l p_k|i p_k|j ((d_kl^2)/2)}, where \emph{S} is  the  number  of  species, \eqn{p_{k|i}}{p_k|i} is the relative abundance of species \emph{k} within community \emph{i}; \eqn{d_{kl}}{d_kl} is the (phylogenetic or functional) dissimilarity between species \emph{k} and \emph{l}.}
  \item{scannf}{
a logical value indicating whether the screeplot (with eigenvalues) should be displayed.}
  \item{nf}{
if \code{scannf} is \code{FALSE}, an integer indicating the number of kept axes.}
  \item{w}{
either a string or a numeric vector of positive values that indicates how the rows of \code{df} (the communities) should be weighted. If \code{w="classic"}, the weights are defined from the sum of the values in each row (e.g. sum of all species abundances within a community). If \code{w="independence"}, then the weight attributed to a row of \code{df} (a community) is the product of the weight attributed to a level of factor A with the weight attributed to a level of factor B. If a vector of strings is given, only the first one is retained. If numeric, values in \code{w} must be in the same order as the rows of \code{df} (see Pavoine et al. 2013 for details on the definition of these weights).}
\item{tol}{
a numeric tolerance threshold: a value between -\code{tol} and \code{tol} is considered as null.}
}
\value{
The functions \code{crossdpcoa_maineffect}, \code{crossdpcoa_version1} and \code{crossdpcoa_version2} return a list containing the following information used for computing the crossed-DPCoA:
\item{l1}{coordinates of the columns of \code{df} (the species).}
\item{l2}{coordinates of the levels of factor A.}
\item{l3}{(for functions \code{crossdpcoa_version1} and \code{crossdpcoa_version2} only) coordinates of the rows of \code{df} (the communities).}
\item{eig}{the eigenvalues.}
\item{lX}{the weights attributed to the columns of \code{df} (species).}
\item{lA}{the weights attributed to the levels of factor A.}
\item{lB}{the weights attributed to the levels of factor B.}
\item{lC}{the weights attributed to the rows of \code{df} (communities).}
\item{div}{a numeric vector with the apportionment of Rao's quadratic diversity (APQE).}
\item{call}{the \code{call} function.}
}
\references{
Pavoine, S., Blondel, J., Dufour, A.-B., Gasc, A., Bonsall, M.B. (2013) A new technique for analysing interacting factors affecting biodiversity patterns: crossed-DPCoA. \emph{PloS One}, \bold{8}, e54530. 
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
\dontrun{
if(require(ape) && require(phylobase) && require(adephylo) 
   && require(adegraphics)){
O <- adegpar()$plabels$optim
adegpar("plabels.optim" = TRUE)

data(birdData)
phy <- read.tree(text=birdData$tre)
phydis <- sqrt(distTips(phy, method="nNodes")+1)

fau <- birdData$fau[, phy$tip.label]
facA <- birdData$facA
facB <- birdData$facB

#Here factor B is put first to analyze 
#the main effect of the strata:
cd_mainB <- crossdpcoa_maineffect(fau, facB, facA, phydis, w=rep(1/30, 30), scannf = FALSE)
barplot(cd_mainB$eig)
cd_mainB$eig[1:2]/sum(cd_mainB$eig)

#Positions of the levels of factor B on its principal axes:
s.label(cd_mainB$l2)
# The "d" value on graphs indicates the length of the edge of a grid cell (scale of the graphic). 

#The coordinates of the species on the same axes 
# can be displayed in front of the phylogeny 
# (several possibilities are provided below, 
# the last one use package adephylo)):
mainBl1.4d <- phylo4d(phy, as.matrix(cd_mainB$l1))
dotp4d(mainBl1.4d, center = FALSE, scale = FALSE)
barp4d(mainBl1.4d, center = FALSE, scale = FALSE)
gridp4d(mainBl1.4d, center = FALSE, scale = FALSE)
parmar <- par()$mar
par(mar=rep(.1,4))
table.phylo4d(mainBl1.4d, show.node=FALSE, symbol="squares",
    center=FALSE, scale=FALSE, cex.label=0.5, ratio.tree=0.7)
par(mar=parmar)

#If factor A is put first, the analysis focus 
#on the main effect of the region:
cd_mainA <- crossdpcoa_maineffect(fau, facA, facB, phydis, w=rep(1/30, 30), scannf = FALSE)
barplot(cd_mainA$eig)
cd_mainA$eig[1:2]/sum(cd_mainA$eig)

#Positions of the levels of factor A on its principal axes:
s.label(cd_mainA$l2)
# The "d" value on graphs indicates the length of the edge of a grid cell (scale of the graphic). 

#The coordinates of the species on the same axes 
# can be displayed in front of the phylogeny
# (several possibilities are provided below, 
# the last one use package adephylo)):
mainAl1.4d <- phylo4d(phy, as.matrix(cd_mainA$l1))
dotp4d(mainAl1.4d, center = FALSE, scale = FALSE)
barp4d(mainAl1.4d, center = FALSE, scale = FALSE)
gridp4d(mainAl1.4d, center = FALSE, scale = FALSE)
parmar <- par()$mar
par(mar=rep(.1,4))
table.phylo4d(mainAl1.4d, show.node=FALSE, symbol="squares", 
   center=FALSE, scale=FALSE, cex.label=0.5, ratio.tree=0.7)
par(mar=parmar)

#Crossed DPCoA Version 1
cd_v1 <- crossdpcoa_version1(fau, facA, facB, phydis, w=rep(1/30, 30), scannf = FALSE)
#Proportion of SS(A) expressed by the two first axes:
cd_v1$eig[1:2]/sum(cd_v1$eig)
#To view the positions of the locations on the first two axes, write:
s.label(cd_v1$l2)
#To view the positions of all communities on the first two axes, write:
s.label(cd_v1$l3)
#To view the positions of the species on the first two axes in front of the phylogeny, write:
v1l1.4d <- phylo4d(phy, as.matrix(cd_v1$l1))
# (then several functions can be used as shown below, 
# the last function, table.phylo4d, is from package adephylo)):
dotp4d(v1l1.4d, center = FALSE, scale = FALSE)
barp4d(v1l1.4d, center = FALSE, scale = FALSE)
gridp4d(v1l1.4d, center = FALSE, scale = FALSE)
parmar <- par()$mar
par(mar=rep(.1,4))
table.phylo4d(v1l1.4d, show.node=FALSE, symbol="squares", 
   center=FALSE, scale=FALSE, cex.label=0.5, ratio.tree=0.7)
par(mar=parmar)

#Crossed DPCoA Version 2
#Crossed DPCoA version 2 can now be performed as follows:
cd_v2 <- crossdpcoa_version2(fau, facA, facB, phydis, w=rep(1/30, 30), scannf = FALSE)
#Proportion of variation among levels of factor A 
#in the subspace orthogonal to the principal axes of B
#expressed by the two first axes:
cd_v2$eig[1:2]/sum(cd_v2$eig)
#To view the positions of the locations on the first two axes, write:
s.label(cd_v2$l2)
#To view the positions of all communities on the first two axes, write:
s.label(cd_v2$l3)
#To view the positions of the species on the first two axes in front of the phylogeny, write:
v2l1.4d <- phylo4d(phy, as.matrix(cd_v2$l1))
# (then several functions can be used as shown below, 
# the last function, table.phylo4d, is from package adephylo)):
dotp4d(v2l1.4d, center = FALSE, scale = FALSE)
barp4d(v2l1.4d, center = FALSE, scale = FALSE)
gridp4d(v2l1.4d, center = FALSE, scale = FALSE)
parmar <- par()$mar
par(mar=rep(.1,4))
table.phylo4d(v2l1.4d, show.node=FALSE, symbol="squares", 
   center=FALSE, scale=FALSE, cex.label=0.5, ratio.tree=0.7)
par(mar=parmar)

adegpar("plabels.optim" = O)
}
}
}
\keyword{models}
