\name{batcomm}
\alias{batcomm}
\docType{data}
\title{
Bat Abundance and Phylogeny Along a Disturbance Gradient in a Neotropical Rainforest
}
\description{
The data were collected by Medellin et al. (2000) on bats in four habitats in the Selva Lacandona of Chiapas, Mexico. The phylogeny is Fritz et al. (2009) phylogeny pruned for retaining only the species present in Medellin et al. data set.
}
\usage{data("batcomm")}
\format{
\code{batcomm} is a list of two components:

\code{ab}, a data frame with habitats as rows, species as columns and abundance of species in habitats as entries.

\code{tre}, a phylogenetic tree in newick format for all bat species in \code{ab}.
}
\details{
In \code{ab}, four habitat types were analyzed: F=rainforest; P=cacao plantations; O=oldfields; C=cornfields (Pavoine 2016).

The resolution of the bat phylogeny is uncertain especially at the older nodes (Pavoine 2016).
}
\source{
\url{http://www.oikosjournal.org/appendix/oik-03262}
}
\references{
Medellin, R., Equihua M., Amin,  M.A. (2000) Bat diversity and abundance as indicators of disturbance in Neotropical rainforest. \emph{Conservation Biology}, \bold{14}, 1666--1675.

Fritz, S.A., Bininda-Emonds, O.R.P., Purvis, A. (2009) Geographic variation in predictors of mammalian extinction risk: big is bad, but only in the tropics. \emph{Ecology Letters}, \bold{12}, 538--549.

Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 
}
\examples{
\dontrun{
if(require(ape)){
data(batcomm)
phy <- read.tree(text=batcomm$tre)
plot(phy)
ab <- batcomm$ab[,phy$tip.label]
plot(abgevodivparam(phy, ab, q=0:4))
}
}
}
\keyword{datasets}
