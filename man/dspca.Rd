\name{dspca}
\alias{dspca}
\alias{plot.dspca}
\title{
Double similarity principal component analysis
}
\description{
The function \code{dspca} performs the double similarity principal component analysis (DSPCA) (Pavoine 2019): an ordination approach to analyse functional or phylogenetic similarities between species communities. The function \code{plot.dspca} displays factorial maps of dspca.
}
\usage{
dspca(comm, S=NULL, tol=1e-8)

\method{plot}{dspca}(x, xaxis = 1, yaxis = 2, labels = TRUE,
    arrows = TRUE, points = FALSE, autolab = TRUE, title = NULL,
    colors = NULL, type = c("X&Y","X","Y"), zoom = TRUE, \dots)
}
\arguments{
  \item{comm}{
a data frame or a matrix typically with communities (or sites, plots, etc.) as rows, species as columns and presence/absence (1/0) or an index of abundance as entries. 
}
  \item{S}{
a matrix of functional or phylogenetic similarities between species (species as rows and columns in the same order as in comm).
}
  \item{tol}{
a numeric tolerance threshold: a value between -tol and tol is considered as zero.
}

\item{x}{
an object of class \code{dspca}.
}
  \item{xaxis}{
the number of the DSPCA axis chosen for the x-axis of the 3d plot.
}
  \item{yaxis}{
the number of the DSPCA axis chosen for the y-axis of the 3d plot.
}
\item{labels}{
a logical specifying whether the names of the species and those of the communities must be displayed in the factorial maps.
}
\item{arrows}{
a logical specifying whether arrows must be used to indicate the positions of the species and those of the communities.
}
\item{points}{
a logical specifying whether points must be used to indicate the positions of the species and those of the communities.
} 
\item{autolab}{
a logical indicating if a function autoLab extracted from package \code{FactoMineR} version 2.6 must be used to display the labels. If there are more than 50 species or more than 50 communities, this option may be time-consuming. If \code{autolab = FALSE}, the names of species and those of the communities will be placed quickly but may overlap.
}
\item{title}{
a string (if \code{type="X"} or \code{type="Y"}) or a vector of strings (if \code{type="X&Y"}) giving the title of each graph.
}
\item{colors}{
a vector (if \code{type="X"} or \code{type="Y"}) providing the color or the colors to be used for species names, points and arrows (if \code{type="X"}), or community names, points and arrows (if \code{type="Y"}), OR, if \code{type="X&Y"}, a list of two vectors, the first one with the colors to be used for species names, points and arrows and the other one for community names, points and arrows.  
}
\item{type}{
a string that can be \code{"X&Y"},\code{"X"},\code{"Y"}. If \code{type="X"}, species scores are displayed. If \code{type="Y"}, community scores are displayed. If \code{"X&Y"}, both species and community scores are displayed.
}
\item{zoom}{
a logical indicating if the graphs must zoom on the positions of the species and those of the communities (defaults to \code{TRUE}). If \code{zoom=FALSE} then species and communities are positioned in the full circle of unit radius.
}
\item{\dots}{
further arguments passed to or from other methods. These must be common to functions plot, arrows, points, and autoLab or text.
}

}
\value{
\code{dspca} returns a list of the following objects:
\item{eig}{Final eigenvalues: positive eigenvalues of the matrix of similarities among communities.}
\item{X}{final coordinates of the species (\eqn{\mathbf{X}_{final}}{X_final}): matrix with the coordinates of the species on the principal components associated with the matrix of similarities among communities. The names of the matrix start with "CPC" indicating "communities' principal component".}
\item{Y}{final coordinates of the communities (\eqn{\mathbf{Y}_{final}}{Y_final}): matrix with the coordinates of
the communities on the principal components associated with the matrix of similarities among communities. The names of the matrix start with "CPC" indicating "communities' principal component".} 
\item{Scom}{matrix of similarities among communities (obtained with coefficient SOchiai (Pavoine and Ricotta 2014))}.
}
\details{
Coordinates can be visualized with graphic tools available in R. Examples are provided below (see section named "examples").
}
\references{
Pavoine, S. (2019) An ordination approach to explore similarities among communities. \emph{Journal of Theoretical Biology}, \bold{462}, 85--96.

Pavoine, S. and Ricotta, C. (2014) Functional and phylogenetic similarity among communities. \emph{Methods in Ecology and Evolution}, \bold{5}, 666--675.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evoNSCA}}, \code{\link{evopcachord}}, \code{\link{evopcahellinger}}, \code{\link{dsimcom}}
}
\examples{
\dontrun{
if(require(ade4)){

data(macroloire, package="ade4") 
# Analysis of the feeding habits of macroinvertebrates in the Loire river, France

# The matrix named feed below contains feeding attributes as rows,
# species as columns, and affinities (proportions) as entries.
feed <- macroloire$traits[, -(1:4)]

# Feeding habits comprise seven categories: engulfers, shredders, scrapers,
# deposit-feeders, active filter-feeders, passive filter-feeders and piercers, in this order.

# Functional similarities among species are computed as indicated in Pavoine and Ricotta (2014)
s_species_ochiai_feed <- dsimFun(feed, vartype = "P", method=4)

# To check that this matrix of trait-based similarities 
# among species is positive semidefinite (p.s.d.)
# we have to verify that its eigenvalues are all nonnegative:
all(eigen(s_species_ochiai_feed)$val>-(1e-8))

dspca_feed <- dspca(t(macroloire$fau), s_species_ochiai_feed)

# eigenvalue distribution:
barplot(dspca_feed$eig)
# The eigenvalues show strong similarities in the composition of the sites 
# regarding the feeding habits of the macroinvertebrates.

# In this data set, communities (sites) are named from S1 to S38 from upstream to downstream
rownames(t(macroloire$fau))

# species are coded from E1 to E40
colnames(t(macroloire$fau))

plot(dspca_feed,  autolab=TRUE, 
    colors= list(ifelse(feed[,6]>0.5,"red","black"), 
    hcl.colors(38, "cyan")))

# The graphical display shows that most sites are characterized by the dominance of a few species 
# (which are mostly passive filter-feeders, in red). 
# In a few sites mostly upstream (in blue), the diversity in feeding habits is higher. 
}
}
}
\keyword{models}




