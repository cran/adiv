\name{birdData}
\alias{birdData}
\docType{data}
\title{
Avian Communities along Successional Forest Gradients
}
\description{
The density of bird species in 3 locations in the Mediterranean region and 2 locations in the central Europe region were collected by Blondel and Farre (1988) along a gradient of habitats.
The densities and a phylogenetic tree are provided for all observed species.
}
\usage{data("birdData")}
\format{
\code{birdData} is a list of 4 objects:

\code{fau}, a data frame with communities as rows, species as columns and densities as entries. The names of the communities start with the first three letters of the location and end with the seral stage number. See details.

\code{tre}, a string. It contains the phylogeny in a newick format.

\code{facA}, a factor. It indicates which location each community belongs to.

\code{facB}, a factor. It indicates which seral stage each community belongs to.
}
\details{
The locations are Alg = Algeria, Bur = Burgundy, Cor = Corsica Island, Pol = Poland, Pro = Provence.

Three are in the Mediterranean region (Provence, southern France; Corsica Island, southern France; and north east Algeria) and two in the central European region (Burgundy, central France; and Poland). 

In each location, a habitat gradient has been conventionally divided into six seral stages (intermediate stages found in forest ecosystems advancing towards their climax stage after a disturbance event) in such a way that all five selected habitat gradients match one another reasonably well in terms of the number, patterns and overall structure of habitats: from low bushy vegetation, less than 1 m height (stage 1), to forests with trees at least 20 m high (stage 6). 

The composite phylogenetic tree was obtained based on Davis' supertree: a strict consensus of 2000 trees (see details in Text S5 of Pavoine et al. 2013).  
}
\source{
Dataset S1 in Pavoine et al. (2013)
}
\references{
Blondel, J. and Farre, H. (1988) The convergent trajectories of bird communities along ecological successions in European forests. \emph{OEcologia (Berlin)}, \bold{75}, 83--93.

Davis, K.E. (2008) \emph{Reweaving the tapestry: a supertree of birds} [PhD thesis]. Glasgow, U.K.: University of Glasgow.

Pavoine, S., Blondel, J., Dufour, A.-B., Gasc, A., Bonsall, M.B. (2013) A new technique for analysing interacting factors affecting biodiversity patterns: crossed-DPCoA. \emph{PloS One}, \bold{8}, e54530. 
}
\examples{
\dontrun{
if(require(ape) && require(adephylo)){
data(birdData)
phy <- read.tree(text=birdData$tre)
phydis <- sqrt(distTips(phy, method="nNodes")+1)

fau <- birdData$fau[, phy$tip.label]
facA <- birdData$facA
facB <- birdData$facB

cd_mainB <- crossdpcoa_maineffect(fau, facB, facA, phydis, w=rep(1/30, 30), scannf = FALSE)
s.label(cd_mainB$l2)
}
}
}
\keyword{datasets}
