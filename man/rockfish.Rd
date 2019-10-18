\name{rockfish}
\alias{rockfish}
\docType{data}
\title{
Rockfish Phylogenetic Diversity in Southern California Bight
}
\description{
The data set was analyzed in Pavoine et al. (2009) to determine whether the observed decrease in the number of individuals and species in the Southern California Bight through years due to fishing activities was accompanied by changes in the phylogenetic structure of the community.
}
\usage{data("rockfish")}
\format{
The format is a list of 3 objects:

\code{tre}, a string: the phylogenetic tree in newick format.

\code{fau}, a data frame with years as rows, species as columns (genus: Sebastes), and abundance as entries.

\code{traits}, a data frame with species as rows (same as in \code{fau}) and two variables as columns: \code{MaxSize}, the maximum body size; and \code{Vul}, an index of species' vulnerability.
}
\details{
The abundance data were obtained from the Marine Recreational Fishery Statistics Survey (MRFSS). We considered the compositions
of rockfish assemblages caught by party and charter boats with hooks and lines from 1980-1986, 1993-1994, 1996, 1998-2007. The phylogenetic tree was obtained from Hyde and Vetter (2007).}
\source{
Appendixes of Pavoine et al. (2009).
}
\references{
Pavoine, S., Love, M., Bonsall, M. (2009) Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: application to rockfish (genus: \emph{Sebastes}) in the Southern California Bight. \emph{Ecology Letters}, \bold{12}, 898--908.

Hyde, J.R. and Vetter, R.D. (2007). The origin, evolution, and
diversification of rockfishes of the genus Sebastes (Cuvier). \emph{Molecular Phylogenetics and Evolution}, \bold{44}, 790--811.
}
\examples{
\dontrun{
if(require(ape)){
data(rockfish)
phy <- read.tree(text=rockfish$tre)
plot(phy)
}
}
}
\keyword{datasets}
