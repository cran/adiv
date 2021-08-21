\name{batcomm}
\alias{batcomm}
\docType{data}
\title{
Bat Abundance and Phylogeny Along a Disturbance Gradient in a Neotropical Rainforest
}
\description{
The data were collected by Medellin et al. (2000) on bats in four habitats in the Selva Lacandona of Chiapas, Mexico. 
Two phylogenetic trees were obtained, one with Fritz et al. (2009) phylogeny pruned for retaining only the species present in Medellin et al. data set; the other one using a consensus tree obtained from 9999 credible trees developed by Upham et al. (2019) and obtained from vertlife.org (see Pavoine and Ricotta (2021) for more details).
}
\usage{data("batcomm")}
\format{
\code{batcomm} is a list of two components:

\code{ab}, a data frame with habitats as rows, species as columns and abundance of species in habitats as entries.

\code{tre}, a phylogenetic tree in newick format for all bat species in \code{ab}.

\code{ab2}, a data frame with habitats as rows, species as columns and abundance of species in habitats as entries.

\code{tre2}, a phylogenetic tree in newick format for all bat species in \code{ab}.
}
\details{
In \code{ab} and \code{ab2}, four habitat types were analyzed: F=rainforest; P=cacao plantations; O=old fields; C=cornfields (Pavoine 2016). 

\code{ab} and \code{tre} contain species names as in Medellin et al. (2000).
\code{ab2} and \code{tre2} contain species names as in Upham et al. (2019).

Synonyms used are: Pteronotus_parnelii / Pteronotus_parnellii; Phyllostomus_stenops / Phylloderma_stenops; Tonatia_brasiliense / Lophostoma_brasiliense; Tonatia_evotis / Lophostoma_evotis; Micronycteris_brachyotis / Lampronycteris_brachyotis; Vampyrodes_major / Vampyrodes_caraccioli. 

\code{tre} was obtained using Fritz et al. (2009) and \code{tre2} using Upham et al. (2019). The resolution of the bat phylogeny obtained from Fritz is uncertain especially at the older nodes (Pavoine 2016).

}
\source{
\url{http://www.oikosjournal.org/appendix/oik-03262}
}
\references{
Fritz, S.A., Bininda-Emonds, O.R.P., Purvis, A. (2009) Geographic variation in predictors of mammalian extinction risk: big is bad, but only in the tropics. \emph{Ecology Letters}, \bold{12}, 538--549.

Medellin, R., Equihua M., Amin,  M.A. (2000) Bat diversity and abundance as indicators of disturbance in Neotropical rainforest. \emph{Conservation Biology}, \bold{14}, 1666--1675.

Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 

Pavoine, S., Ricotta, C. (2021) On the relationships between rarity, uniqueness, distinctiveness, originality and functional/phylogenetic diversity. \emph{Contact authors for information on this study}

Upham, N.S., Esselstyn, J.A., Jetz, W. (2019) Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. \emph{PloS Biology}, \bold{17}, e3000494. 
}
\examples{
\dontrun{
if(require(ape)){
data(batcomm)
phy <- read.tree(text=batcomm$tre)
plot(phy)
ab <- batcomm$ab[,phy$tip.label]
plot(abgevodivparam(phy, ab, q=0:4))

phy2 <- read.tree(text=batcomm$tre2)
plot(phy2)
ab2 <- batcomm$ab2[,phy2$tip.label]
plot(abgevodivparam(phy2, ab2, q=0:4))
}
}
}
\keyword{datasets}
