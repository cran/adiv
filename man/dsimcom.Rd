\name{dsimcom}
\alias{dsimcom}
\alias{SQ}
\alias{dsimTaxo}
\alias{dsimFun}
\alias{dsimTree}
\title{
Functional or Phylogenetic Similarity between Species and Communities
}
\description{
The R function \code{dsimcom} calculates Pavoine and Ricotta (2014) coefficients SSokal-Sneath, SJaccard, SSorensen, SOchiai, and Sbeta of similarities among communities or their complements (1-S, dissimilarities).

Function \code{SQ} calculates Pavoine and Ricotta (2014) index \eqn{S_Q} (similarities between communities) and its complements \eqn{1-S_Q} (dissimilarities between communities).

The functions \code{dsimTax} and \code{dsimTree} calculate pair-wise taxonomic and phylogenetic (dis)similarities between species, respectively. \code{dsimTree} can also be used to calculate pair-wise functional (dis)similarities between species if a functional dendrogram is used to describe species.

Function \code{dsimFun} calculates pair-wise functional (dis)similaties between species.
}
\usage{
dsimcom(comm, Sigma = NULL, method = 1:5, 
option = c("relative", "absolute"), 
type = c("similarity", "dissimilarity"))

SQ(comm, Sigma = NULL, type = c("similarity", "dissimilarity"))

dsimTaxo(tax, method = c(1, 2, 3, 4, 5), 
type = c("similarity", "dissimilarity"))

dsimTree(phyl, method = c(1, 2, 3, 4, 5), rootedge = NULL, 
type = c("similarity", "dissimilarity"))

dsimFun(df, vartype = c("Q", "N", "M", "P"), method = 1:5, 
type = c("similarity", "dissimilarity"))
}
\arguments{
  \item{comm}{a data frame or matrix with communities as rows, species as columns and non-negative values as entries.}
  \item{Sigma}{matrix of similarities among species (species as rows and columns in the same order as in \code{comm}; values in \code{Sigma} are bounded between 0 and 1). The matrix must be nonnegative definite, i.e. all its eigenvalues are nonnegative. Functions \code{dsimTaxo}, \code{dsimTree} and \code{dsimFun} can be used to obtain \code{Sigma} as these functions lead to nonnegative definite matrices.}
  \item{method}{
an  integer  (1, 2, 3, 4, 5)  indicating  which basic coefficient  should  be  used: Sokal-Sneath, Jaccard, Sorensen, Ochiai, beta, respectively.}
  \item{option}{
a  string. If  \code{option = "relative"},  the  rows  of  \code{comm}  are  standardized  into proportions that sum to 1.  If \code{option = "absolute"}, raw values are retained in \code{comm}.}
  \item{tax}{an object of class \code{taxo} (of package ade4).}
  \item{type}{a string. If \code{type = "similarity"}, the functions \code{dsimcom} and \code{SQ} calculates  similarities (S) between communities and the functions \code{dsimTaxo}, \code{dsimTree} and \code{dsimFun} calculates similarities between species. If \code{type = "dissimilarity"}, the functions \code{dsimcom} and \code{SQ} calculates  dissimilarities (1-S) between communities and the functions \code{dsimTaxo}, \code{dsimTree} and \code{dsimFun} calculates dissimilarities between species.}
   \item{phyl}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.}
   \item{rootedge}{a numeric equal to the length of the branch at the nearest common ancestor of  all  species (here referred to as the root node). This branch is thus anterior to the root.  It is ignored if \code{rootedge} is NULL.}
   \item{df}{either a data frame, a matrix or an object of class \code{ktab} of package ade4. If class \code{ktab} is used, each data frame must contain a single type of traits (see argument \code{vartype}).}
   \item{vartype}{a vector of characters indicating the type of traits used in each data frame of \code{df}: \code{"Q"} for quantitative; \code{"N"} for nominal; \code{"M"} for multichoice; \code{"P"} for traits that can be expressed as proportions. Values in type must be in the same order as data frames in \code{df} (one value per table).}
}
\details{
Formulas for the indices are given in Pavoine and Ricotta (2014) main text and appendixes.
}
\value{
If \code{type = "similarities"}, a matrix of similarities between communities.
If \code{type = "dissimilarities"}, an object of class \code{dist} with the dissimilarities between communities.
}
\references{
Pavoine, S., Ricotta, C. (2014) Functional and phylogenetic similarity among communities. \emph{Methods in Ecology and Evolution}, \bold{5}, 666--675.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
\dontrun{
if(require(ade4)){
data(macroloire, package="ade4") 

Ssokalsneath <- dsimcom(t(macroloire$fau), method=1, option=c("relative"))
Sjaccard <- dsimcom(t(macroloire$fau), method=2, option=c("relative"))
Ssorensen <- dsimcom(t(macroloire$fau), method=3, option=c("relative"))
Sochiai <- dsimcom(t(macroloire$fau), method=4, option=c("relative"))
Sbeta <- dsimcom(t(macroloire$fau), method=5, option=c("relative"))

SQUNIF <- SQ(t(macroloire$fau))

# The taxonomy is contained in macroloire$taxo

s_species_sokalsneath_taxo <- dsimTaxo(macroloire$taxo, method=1) # Using formula a/(a+2b+2c) 
# (see notations below)
s_species_jaccard_taxo <- dsimTaxo(macroloire$taxo, method=2) # Using formula a/(a+b+c) 
s_species_sorensen_taxo <- dsimTaxo(macroloire$taxo, method=3) # Using formula 2a/(2a+b+c) 
s_species_ochiai_taxo <- dsimTaxo(macroloire$taxo, method=4) # Using a/sqrt((a+b)(a+c)) 
s_species_beta_taxo <- dsimTaxo(macroloire$taxo, method=5) # Using 4a/(4a+b+c)


# To check that these matrices of taxonomic similarities 
# among species are positive semidefinite (p.s.d.)
# we have to verify that their eigenvalues are all nonnegative:
all(eigen(s_species_sokalsneath_taxo)$val>-(1e-8))
all(eigen(s_species_jaccard_taxo)$val>-(1e-8))
all(eigen(s_species_sorensen_taxo)$val>-(1e-8))
all(eigen(s_species_ochiai_taxo)$val>-(1e-8))
all(eigen(s_species_beta_taxo)$val>-(1e-8))

# Compositions of the communities
m <- t(macroloire$fau[rownames(s_species_sokalsneath_taxo), ])
# Taxomonic similarities among species
s <- s_species_sokalsneath_taxo

Ssokalsneath_taxo <- dsimcom(m, s, method = 1)
SQwith_species_sokalsneath_taxo <- SQ(m, s)

# Compositions of the communities
m <- t(macroloire$fau[rownames(s_species_jaccard_taxo), ])
# Taxonomic similarities among species
s <- s_species_jaccard_taxo  
Sjaccard_taxo <- dsimcom(m, s, method = 2)
SQwith_species_jaccard_taxo <- SQ(m, s)

# Compositions of the communities
m <- t(macroloire$fau[rownames(s_species_sorensen_taxo), ])
# Taxonomic similarities among species
s <- s_species_sorensen_taxo  

Ssorensen_taxo <- dsimcom(m, s, method = 3)
SQwith_species_sorensen_taxo <- SQ(m, s)

# Compositions of the communities
m <- t(macroloire$fau[rownames(s_species_ochiai_taxo), ])
# Taxonomic similarities among species
s <- s_species_ochiai_taxo  

Sochiai_taxo <- dsimcom(m, s, method = 4)
SQwith_species_ochiai_taxo <- SQ(m, s)


# Compositions of the communities
m <- t(macroloire$fau[rownames(s_species_beta_taxo), ])
# Taxonomic similarities among species
s <- s_species_beta_taxo  

Sbeta_taxo <- dsimcom(m, s, method = 5)
SQwith_species_beta_taxo <- SQ(m, s)

# The matrix named feed below contains feeding attributes as rows,
# species as columns, and affinities (proportions) as entries.
feed <- macroloire$traits[ ,-(1:4)]

# Feeding habits comprise seven categories: engulfers, shredders, scrapers,
# deposit-feeders, active filter-feeders, passive filter-feeders and piercers, in this order.

# Functional similarities among species are computed as indicated in the main text
s_species_sokalsneath_feed <- dsimFun(feed, vartype = "P", method=1)
s_species_jaccard_feed <- dsimFun(feed, vartype = "P", method=2)
s_species_sorensen_feed <- dsimFun(feed, vartype = "P", method=3)
s_species_ochiai_feed <- dsimFun(feed, vartype = "P", method=4)
s_species_beta_feed <- dsimFun(feed, vartype = "P", method=5)

all(eigen(s_species_sokalsneath_feed)$val>-(1e-8))
all(eigen(s_species_jaccard_feed)$val>-(1e-8))
all(eigen(s_species_sorensen_feed)$val>-(1e-8))
all(eigen(s_species_ochiai_feed)$val>-(1e-8))
all(eigen(s_species_beta_feed)$val>-(1e-8))

Ssokalsneath_feed <- dsimcom(t(macroloire$fau), s_species_sokalsneath_feed, method=1)
SQwith_species_sokalsneath_feed <- SQ(t(macroloire$fau), s_species_sokalsneath_feed)

Sjaccard_feed <- dsimcom(t(macroloire$fau), s_species_jaccard_feed, method=2)
SQwith_species_jaccard_feed <- SQ(t(macroloire$fau), s_species_jaccard_feed)

Ssorensen_feed <- dsimcom(t(macroloire$fau), s_species_sorensen_feed, method=3)
SQwith_species_sorensen_feed <- SQ(t(macroloire$fau), s_species_sorensen_feed)

Sochiai_feed <- dsimcom(t(macroloire$fau), s_species_ochiai_feed, method=4)
SQwith_species_ochiai_feed <- SQ(t(macroloire$fau), s_species_ochiai_feed)

Sbeta_feed <- dsimcom(t(macroloire$fau), s_species_sorensen_feed, method=5)


all(eigen(Ssokalsneath_feed)$val>-(1e-8))
all(eigen(Sjaccard_feed)$val>-(1e-8))
all(eigen(Ssorensen_feed)$val>-(1e-8))
all(eigen(Sochiai_feed)$val>-(1e-8))
all(eigen(Sbeta_feed)$val>-(1e-8))


par(mfrow=c(3, 5))
plot(SQUNIF, Ssokalsneath, xlim=c(0,1), ylim=c(0,1), asp=1)
segments(0, 0, 1,1)
plot(SQUNIF, Sjaccard, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQUNIF, Ssorensen, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQUNIF, Sochiai, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQUNIF, Sbeta, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)

plot(SQwith_species_sokalsneath_taxo, Ssokalsneath_taxo, xlim=c(0,1), ylim=c(0,1), asp=1)
segments(0, 0, 1,1)
plot(SQwith_species_jaccard_taxo, Sjaccard_taxo, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQwith_species_sorensen_taxo, Ssorensen_taxo, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQwith_species_ochiai_taxo, Sochiai_taxo, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQwith_species_beta_taxo, Sbeta_taxo, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)

plot(SQwith_species_sokalsneath_feed, Ssokalsneath_feed, xlim=c(0,1), ylim=c(0,1), asp=1)
segments(0, 0, 1,1)
plot(SQwith_species_jaccard_feed, Sjaccard_feed, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQwith_species_sorensen_feed, Ssorensen_feed, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQwith_species_ochiai_feed, Sochiai_feed, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)
plot(SQwith_species_sorensen_feed, Sbeta_feed, xlim=c(0,1), ylim=c(0,1) , asp=1)
segments(0, 0, 1,1)

par(mfrow=c(1,1))
}
}
}
\keyword{models}
