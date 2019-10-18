\name{aptree}
\alias{apd}
\alias{aptree}
\alias{plot.aptree}
\alias{abgaptree}
\alias{rtestaptree}
\alias{plot.rtestaptree}
\alias{tecAptree}
\alias{pIa}
\title{
Hierarchical Partitioning of Evolutionary and Ecological Patterns in the Organization of Phylogenetically-Structured Species Assemblages
}
\description{
\code{apd} performs Hardy (2008)'s test for phylogenetic structure in species abundance distribution; 

\code{aptree} apportions the diversity (according to index \eqn{I_a}{Ia} by Pavoine et al. 2009) within one or several communities between evolutionary periods; 

\code{plot.aptree} displays the phylogenetic tree with vertical lines at each speciation event (limits of the evolutionary periods), the first period starts at the tips and the last one ends at the root node; the phylogenetic tree is pruned retaining only the species present in at least one of the observed communities;

\code{abgaptree} provides the apportionment of alpha, beta and gamma diversities between evolutionary periods, according to index \eqn{I_a}{Ia} by Pavoine et al. (2009);

\code{rtestaptree} performs the test of phylogenetic signal in the differences between communities at each evolutionary periods; 

\code{plot.rtestaptree} displays the phylogenetic tree with vertical lines at each speciation event (limits of the evolutionary periods; see above); colours are used to highlight the periods where the dissimilarities between communities are different from that expected at random;
 
\code{tecAptree} provides technical information for the apportionment of diversity between evolutionary periods; 

\code{pIa} calculates the index \eqn{I_a}{Ia} by Pavoine et al. (2009) within each community.
}
\usage{
apd(phyl, comm, wcom = c("even", "speciesab"), nrep = 99, 
alter = "two-sided", tol = 1e-08)

aptree(phyl, comm, exponent = 2, tol = 1e-08)

\method{plot}{aptree}(x, col.line = 'blue', ...)

abgaptree(phyl, comm, exponent = 2, 
wcom = c("even", "speciesab"), tol = 1e-08)

rtestaptree(phyl, comm, nrep = 99, alter = "two-sided", 
exponent = 2, wcom = c("even", "speciesab"), tol = 1e-08)

\method{plot}{rtestaptree}(x, col.line = c("blue", "red"), 
alpha = 0.05, ...)

tecAptree(phyl, v = NULL, tol = 1e-08)

pIa(phyl, comm, exponent = 2, tol = 1e-08)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.
}
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and presence/absence or an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips. In function \code{aptree}, \code{comm} can be a vector that provides the presence/absence or an index of abundance within a single community.
}
  \item{wcom}{
a numeric vector that gives the weight attributed to the community. The weights must be positive and their sum equals 1.
}
  \item{nrep}{
a numeric that gives the number of permutations.
}
  \item{alter}{
a string specifying the alternative hypothesis; it must be one of "greater", "less" or "two-sided".
}
  \item{tol}{
a numeric. If the absolute value of a statistic is less than tol, this statistic is considered equal to zero.
}
  \item{exponent}{
a numeric that gives the value of parameter \emph{a} in index \eqn{I_a}{Ia}.
}
  \item{x}{
in \code{plot.aptree}, \code{x} is an object inheriting class \code{aptree} obtained with function \code{aptree}. In \code{plot.rtestaptree}, \code{x} is an object inheriting class \code{rtestaptree} obtained with function \code{rtestaptree}.
}
  \item{col.line}{
in \code{plot.aptree}, \code{col.line} is a string which attributes a colour to the vertical lines placed at each speciation event and defining the evolutionary periods. In \code{plot.rtestaptree}, it is a vector with two strings. These strings give colours to the vertical lines placed at each speciation event. The first colour is used when the differences between communities at the evolutionary period are not significantly different from random; the second colour is used when they are. 
}
  \item{alpha}{
a numeric: the nominal alpha level for significancy (the p-values calculated with function \code{rtestaptree} are compared to alpha to determine the result of the test).
}
  \item{\dots}{
further arguments passed to or from other methods.
}
  \item{v}{
either \code{NULL} or a vector that provides the presence/absence or an index of abundance of species within a single community.
}
}
\details{
The approaches developed in these functions rely on a parametric index of phylogenetic diversity named \eqn{I_a}{Ia}. The parameter \emph{a} controls the importance given to rare versus abundant species in communities. Index \eqn{I_a}{Ia} generalizes Rao's quadratic entropy (QE) applied to phylogenetic distances between species (when \emph{a}=2) and Faith's Phylogenetic Diversity index (PD) (when \emph{a}=0). When \emph{a} tends towards 1, the index is a generalization of the Shannon index of diversity applied to phylogenetic data in addition to abundance data. In Pavoine et al. (2009), we developed this index and demonstrated how it can be used to partition diversity simultaneously across evolutionary periods in the phylogeny and across spatial (e.g. local communities in a region) and/or time units (e.g. a community investigated yearly).
}
\value{
The function \code{apd} returns an object of class \code{randtest} with the results of the test (see function \code{randtest} in package ade4).

The function \code{aptree} returns a data frame with the evolutionary periods as rows, the communities as columns and the diversity values as entries.

The function \code{plot.aptree} returns a graph.

The function \code{abgaptree} returns a data frame with the evolutionary periods as rows, alpha diversity, beta diversity and gamma diversity as columns and the diversity values as entries.

The function \code{rtestaptree} returns an object of class \code{krandtest} with the results of the permutation tests. (see function \code{krandtest} in package ade4)

The function \code{plot.rtestaptree} returns a graph.

The function \code{tecAptree} returns a list. If \code{v} is \code{NULL}, the values of the list are:
\item{h}{the height at which each evolutionary period ends;}
\item{plength}{period length;} 
\item{ngroups}{number of monophyletic groups per evolutionary period;}
\item{list}{list of the species per monophyletic group at each evolutionary period;}
\item{call}{original call.}
If \code{v} contains a vector of presence/absence or abundance, the following object is added in the output: 
\item{relab}{the relative abundance (sum of species' presences or abundances depending on \code{v}) of each monophyletic group at each evolutionary period.}

The function \code{pIa} returns a data frame with communities as rows and only one column. This column gives, for each community, the value taken by index \eqn{I_a}{Ia} of phylogenetic diversity developed by Pavoine et al. (2009).  
}
\references{
Pavoine, S., Love, M., Bonsall, M.B. (2009) Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: application to rockfish (genus: Sebastes) in the Southern California Bight. \emph{Ecology Letters}, \bold{12}, 898--908.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
with contributions of Stephane Dray.
}
\seealso{
\code{\link{QE}}
}
\examples{
\dontrun{
if(require(ape)){
data(rockfish)
phy <- read.tree(text=rockfish$tre)
ABG <- abgaptree(phy, rockfish$fau, wcom="speciesab")
colSums(ABG)
A <- aptree(phy, rockfish$fau)
colSums(A)
plot(A, cex=0.5)
P <- pIa(phy, rockfish$fau)
P
T <- apd(phy, rockfish$fau)
plot(T)
#R <- rtestaptree(phy, rockfish$fau, nrep=999, wcom="speciesab")
#plot(R)
TA <- tecAptree(phy)
TA$h
}
}
}
\keyword{models}
