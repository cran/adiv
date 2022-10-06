\name{rlqESLTP}
\alias{rlqESLTP}
\alias{plot.rlqESLTP}
\alias{summary.rlqESLTP}
\title{
Linking Patterns in Phylogeny, Traits, Abiotic Variables
and Space}
\description{
An extention of the RLQ approach to identify potential environmental filters in species' traits in a phylogenetic and spatial context. 
}
\usage{
rlqESLTP(dudiE, dudiS, dudiL, dudiT, dudiP, ...)

\method{plot}{rlqESLTP}(x, which = NULL, phyl = NULL, 
    xy = NULL, traits = NULL, env = NULL, type = NULL, 
    ax = 1, disp = c("dots", "bars", "grid"), ...)

\method{summary}{rlqESLTP}(object, ...)
}
\arguments{
  \item{dudiE}{
an object of class \code{dudi} of package ade4 (e.g. a principal component applied to environmental variables or a principal coordinate analysis applied to environmental distances between plots). Plots must be weighted as in \code{dudiL}. \code{dudiE} can be set to \code{NULL} if \code{dudiS} is not.
}
  \item{dudiS}{
an object of class \code{dudi} of package ade4 (e.g. an eigenvector decomposition applied to a spatial neighbourhood matrix (see Pavoine et al. 2011 for exemples)). Plots must be weighted as in \code{dudiL}. \code{dudiS} can be set to \code{NULL} if \code{dudiE} is not.
}
  \item{dudiL}{
an object of class \code{dudi} obtained with function \code{dudi.coa} of ade4 applied to a plot (rows) by species (columns) data frame. \code{dudiL} cannot be set to \code{NULL}.
}
  \item{dudiT}{
an object of class \code{dudi} of package ade4 (e.g. a principal component analysis applied to functional traits of species or a principal coordinate analysis applied to functional distances between species). Species must be weighted as in \code{dudiL}. \code{dudiT} can be set to \code{NULL} if \code{dudiP} is not.
}
  \item{dudiP}{
an object of class \code{dudi} of package ade4 (e.g. a principal coordinate analysis applied to phylogenetic distances between species). Species must be weighted as in \code{dudiL}. \code{dudiP} can be set to \code{NULL} if \code{dudiT} is not.
}
  \item{x}{
an object of class \code{rlqESTP} (obtained with function \code{rqlESTP}).}
  \item{object}{
an object of class \code{rlqESTP} (obtained with function \code{rqlESTP}).}
  \item{which}{
a character that might be \code{"E"} (in that case arguments \code{env} and \code{type} must be given); \code{"S"} (in that case argument \code{xy} must be given); \code{"T"} (in that case arguments \code{traits} and \code{type} must be given); or \code{"P"} (in that case argument \code{phyl} must be given).}
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}.
}
  \item{xy}{
a data frame with two columns giving the coordinates of the plots (longitude, latitude in that order).}
  \item{traits}{
a list of data frames for the traits, each data frame contains variables from a single statistical type (see argument \code{type} below). NAs are allowed. However, in that case, the calculations depend, for each trait, on the subset of available data only.}
  \item{env}{
a list of data frames for the environmental variables, each data frame contains variables from a single statistical type (see argument \code{type} below). NAs are allowed. However, in that case, the calculations depend, for each variable, on the subset of available data only.}
  \item{type}{
a vector that provides the type of each table to be analysed in \code{env} (if environemental variables should be displayed) or \code{traits} (if the focus is on species' traits). The possible types are \code{"Q"} (quantitative, numeric), \code{"O"} (ordinal), \code{"N"} (nominal), \code{"D"} (dichotomous), \code{"F"} (fuzzy, or expressed as a proportion), \code{"B"} (multichoice nominal variables, coded by binary columns), \code{"C"} (circular). Values in type must be in the same order as in \code{env} or \code{traits}.}
  \item{ax}{
a numeric indicating the axis of interest.}
  \item{disp}{
a string indicating which method to use to display species scores in front of the phylogeny: one of \code{"dots"} (for a dot plot), \code{"bars"} (for a barplot), and \code{"grid"} (for a gridplot for vizualizing scores by a color gradient).}
  \item{\dots}{
further arguments passed to or from other methods.}
}
\value{
Function \code{rlqESLTP} returns an object of classes \code{rlqESLTP}, \code{rlq} and \code{dudi}.
It is a list of 26 objects:
\item{tab}{a data frame. Crossed Table (CT): crossing the columns of the merged trait and phylogenetic table with those of the merged environmental and spatial table.}
\item{cw}{a vector of numerics. Weights attributed to the columns of the merged trait and phylogenetic table.}
\item{lw}{a vector of numerics. weights attributed to the columns of the merged environmental and spatial table.}
\item{eig}{a vector of numerics. The vector of eigenvalues.}
\item{rank}{an integer. The total number of axes in the analysis.}
\item{nf}{a numeric. The number of axes kept.}
\item{c1}{a data frame. Principal axes. Normed scores for the columns of merged trait and phylogenetic table.}
\item{co}{a data frame. Scores for the columns of merged trait and phylogenetic table.}
\item{l1}{a data frame. Principal axes. Normed scores for the columns of merged environmental and spatial table.}
\item{li}{a data frame. Scores for the columns of merged environmental and spatial table.}
\item{call}{call}
\item{lQ}{data frame. Scores for the species.}       
\item{mQ}{data frame. Normed scores for the species.}
\item{lR}{data frame. Scores of the plots.}
\item{mR}{data frame. Normed scores for the plots.}
\item{aQ}{data frame. Correlations trait/phylogeny axes / coinertia axes.}
\item{aR}{data frame. Correlations environment/space axes / coinertia axes.}
\item{lR_givenE}{data frame. Contributions of environmental information to plot scores.}
\item{lR_givenS}{data frame. Contributions of spatial information to plot scores.}
\item{lQ_givenT}{data frame. Contributions of trait information to species scores.}
\item{lQ_givenP}{data frame. Contributions of phylogenetic information to species scores.}
\item{row.w}{a vector. Weights attributed to plots.}
\item{col.w}{a vector. Weights attributed to species.}
\item{dudiL}{object of class \code{dudi} (\code{dudiL}).}    
\item{dudiR}{object of class \code{dudi} (merging information on environment and space).}
\item{dudiQ}{object of class \code{dudi} (merging information on traits and phylogeny).} 
}
\references{
Pavoine, S., Vela, E., Gachet, S., de Belair, G., Bonsall, M.B. (2011) Linking patterns in phylogeny, traits, abiotic variables and space: a novel approach to linking environmental filtering and plant community assembly. \emph{Journal of Ecology}, \bold{99}, 165--175.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
\dontrun{
if(require(ade4) && require(adegraphics) && require(ape)){
O <- adegpar()$plabels$optim
adegpar("plabels.optim" = TRUE)

data(mafragh, package="ade4")
xy <- mafragh$xy
#The object that defines the neighbourhood between plots is in
mneig <- mafragh$neig
#mneig is an object of class neig
s.label(mafragh$xy, nb = mafragh$nb, paxes.draw = FALSE)

#The environmental variables are in
env <- mafragh$env
names(env)
#The abundance data are in
flo <- mafragh$flo
# Adjustement of the writing of species names
names <- gsub(" ", "_", mafragh$spenames[colnames(flo), 1])
for(i in 1:26){
names <- gsub(LETTERS[i], letters[i], names)
}
names <- gsub("alisma_plantago", "alisma_plantago_aquatica", names)
colnames(flo) <- names
#The data on traits are in
traits <- lapply(mafragh$traits, function(x) x[colnames(flo), , drop=FALSE])
#The phylogenetic tree is in
phy <- read.tree(text=mafragh$tre)
plot(phy)

#Traits are separated by statistical type. The object traits is a list of data frames.
tabBinary <- prep.binary(traits$tabBinary, c(3, 4))
tabQuantitative <- traits$tabQuantitative
tabCircular <- prep.circular(traits$tabCircular, 1, 12)
tabOrdinal <- traits$tabOrdinal

#Analyses:
coaflo <- dudi.coa(flo, scan = FALSE, nf = 55)
vecspa <- scores.neig(mneig)
pcaspa <- dudi.pca(vecspa, coaflo$lw, scan = FALSE, nf = 96)
#We first removed environmental variables
env <- env[-(8:10)]
env[4:8] <- log(env[4:8])
pcaenv <- dudi.pca(env, coaflo$lw, scale = FALSE, scan = FALSE, nf = 8)
disT <- dist.ktab(ktab.list.df(list(tabBinary,
tabOrdinal[c("Spikiness", "Hairy leaves")])), c("B", "O"))
# The definition of the functional distances between species 
# were slightly different in Pavoine et al. (2011).

pcotraits <- dudi.pco(disT, coaflo$cw, full = TRUE)
pcophy <- dudi.pco(sqrt(as.dist(cophenetic(phy)[names(flo), names(flo)])/2), coaflo$cw, full = TRUE)


rlqmix <- rlqESLTP(pcaenv, pcaspa, coaflo, pcotraits, pcophy, scan = FALSE, nf = 2)
barplot(rlqmix$eig)
rlqmix$eig[1]/sum(rlqmix$eig)
plot(rlqmix, xy=xy, ax=1, wh="S")
plot(rlqmix, phy=phy, ax=1, wh="P")
plot(rlqmix, traits=tabBinary, ax=1, type="B", wh="T")
plot(rlqmix, traits=tabOrdinal[2:3], ax=1, type="O", wh="T")
plot(rlqmix, env=env, ax=1, type="Q", wh="E")

adegpar("plabels.optim" = O)
}
}
}
\keyword{models}
