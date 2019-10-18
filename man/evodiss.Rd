\name{evodiss}
\alias{evodiss}
\title{
Pair-wise Phylogenetic Dissimilarities between Communities
}
\description{
The function calculates PD-dissimilarity indices described and/or discussed in Pavoine (2016). Part of them are parametric indices. 
}
\usage{
evodiss(phyl, comm, method = NULL, q = NULL, 
w = c("evoab", "even", "speciesab"), diag = FALSE, 
upper = FALSE, tol = 1e-08)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}.
}
  \item{comm}{
a data frame typically with communities as rows, species as columns and an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips.
}
  \item{method}{
one of the following string codes: \code{"Minkowski"}, \code{"Euclidean"}, \code{"Manhattan"}, \code{"Chord"}, \code{"ScaledCanberra"}, \code{"Divergence"}, \code{"BC"}, \code{"MH"}, \code{"LG"}, \code{"Hellinger"}, \code{"chi2"}, \code{"Hill"}, \code{"Renyi"}, \code{"C"}, \code{"U"}, \code{"S"}. See Details.
}
  \item{q}{
a vector with nonnegative value(s) for parameter \code{q} only used with \code{method="Minkowski"}, \code{"Hill"}, \code{"Renyi"} (parametric indices). \code{q} controls the importance given to rare compared to abundant evolutionary units. See Details.
}
  \item{w}{
either a numeric vector giving weights for communities (same order as in \code{comm}), or a string among \code{"even"}, \code{"evoab"}, and \code{"speciesab"}. If several strings are given, only the first one is used. \code{w} is used only if \code{method="Hill"}, \code{"Renyi"}, \code{"C"}, \code{"U"}, or \code{"S"} is chosen. See details.
}
  \item{diag}{
logical argument passed to function as.dist (R base).
}
  \item{upper}{
logical argument passed to function as.dist (R base).
}
  \item{tol}{
numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.
}
}
\details{
The indices available are (formulas can be found in Supplementary material Appendix 1 of Pavoine 2016):

\code{"Minkowski"}:\eqn{^qevoD_{Minkowski}}{qevoD_Minkowski}

\code{"Euclidean"}: \eqn{evoD_{Euclidean}}{evoD_Euclidean}

\code{"Manhattan"}: \eqn{evoD_{Manhattan}}{evoD_Manhattan}

\code{Chord}: \eqn{evoD_{Chord}}{evoD_Chord}

\code{ScaledCanberra}: \eqn{evoD_{ScaledCanberra}}{evoD_ScaledCanberra}

\code{Divergence}: \eqn{evoD_{Divergence}}{evoD_Divergence}

\code{BC}: \eqn{evoD_{Bray-Curtis}}{evoD_Bray-Curtis}

\code{MH}: \eqn{evoD_{Morisita-Horn}}{evoD_Morisita-Horn}

\code{LG}: \eqn{evoD_{Profile}}{evoD_Profile}

\code{Hellinger}: \eqn{evoD_{Hellinger}}{evoD_Hellinger}

\code{chi2}: \eqn{evoD_{Chi2}}{evoD_Chi2}

\code{Hill}: \eqn{1-\bar{V}_{q2}}{1-V_q2}

\code{Renyi}: \eqn{^qevoD_{Renyi}}{qevoD_Renyi}

\code{C}: \eqn{1-\bar{C}_{q2}}{1-C_q2}

\code{U}: \eqn{1-\bar{U}_{q2}}{1-U_q2}

\code{S}: \eqn{1-\bar{S}_{q2}}{1-S_q2}

The weights of the communities (argument \code{w}) can be \code{"even"} (even weights, i.e. relative abundances are considered for evolutionary units), \code{"evoab"} (proportional to the summed abundances of all evolutionary units, i.e. absolute abundances are considered for evolutionary units), or \code{"speciesab"} (proportional to the summed abundances of all species). Note that if the phylogenetic tree is ultrametric (the distance from any species to the root is constant), then
options \code{"evoab"} and \code{"speciesab"} are equivalent.}
\value{
An object of class \code{dist} containing the PD-dissimilarities (phylogenetic dissimilarities) between communities.
}
\references{
The methodologies are presented in 

Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 

They gather, in a common framework, indices introduced in  

Chiu, C.-H., Jost, L., Chao, A. (2014) Phylogenetic beta diversity, similarity and differentiation measures based on Hill numbers. \emph{Ecological Monographs}, \bold{84}, 21--44.

and earlier work extended here in a phylogenetic context and reviewed in 

Legendre, P. and De Caceres, M. (2013) Beta diversity as the variance of community data: dissimilarity coefficients and partitioning. \emph{Ecology Letters}, \bold{16}, 951--963.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evodiss_family}}
}
\examples{
\dontrun{
if(require(ape)){

data(batcomm)
phy <- read.tree(text=batcomm$tre)
ab <- batcomm$ab[,phy$tip.label]

evodiss(phy, ab, "LG")
evodiss(phy, ab, "Hellinger")
evodiss(phy, ab, "Chord")

evodiss(phy, ab, "Hill", q=2)
evodiss(phy, ab, "Hill", q=2, w="even")

}
}
}
\keyword{models}
