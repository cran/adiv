\name{evodiv}
\alias{evodiv}
\title{
Indices of Phylogenetic Diversity
}
\description{
The function \code{evodiv} calculates diversity indices that rely on the relative or absolute abundance of features on a phylogenetic tree, with the assumption that the number of features on a given branch of a phylogenetic tree is equal to the length of this branch (see Pavoine 2016).
}
\usage{
evodiv(phyl, comm, method = "full", tol = 1e-8)
}
\arguments{
   \item{phyl}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.}
  \item{comm}{a data frame or a matrix typically with communities as rows, species as columns and an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips. Note that with presence/absence (0/1) data, only feature richness will be calculated correctly.}
  \item{method}{
a string or a vector of strings: one or several of "richness", "GiniSimpson", "Simpson", "Shannon", "Margalef", "Menhinick", "McIntosh", "full". See details.
}
  \item{tol}{
a tolerance threshold (a value between -\code{tol} and \code{tol} is considered equal to zero)}
}
\value{
Function \code{evodiv} returns a matrix with communities as rows and the diversity indices as columns. 
}
\details{
Let \eqn{L_b} the length of branch b in the phylogenetic tree. Let \eqn{a_bi} the absolute abundance of branch \emph{b} in community \emph{i} (sum of abundance of all species descending from it in the phylogenetic tree). \eqn{p_{bi} = a_{bi}/\sum_b L_b a_{bi}}{p_bi = a_bi / (sum_b L_b a_bi)} is the relative abundance of branch \emph{b} in community \emph{i}.
If \code{method="richness"}, the diversity index is the number of features (branch units): \eqn{\sum_b L_b}{sum_b L_b}. It corresponds to Faith (1992) Phylogenetic Diversity index.
If \code{method="GiniSimpson"}, the diversity index is that of Gini (1912) and Simpson (1949): \eqn{1-\sum_b L_b p_{bi}^2}{1-sum_b L_b p_bi^2}.
If \code{method="Simpson"}, the diversity index is that of Simpson (1949): \eqn{1/\sum_b L_b p_{bi}^2}{1/(sum_b L_b p_bi^2)}.
If \code{method="Shannon"}, the diversity index is that of Shannon (1948) with neperian logarithm: \eqn{-\sum_b L_b p_{bi} log(p_{bi})}{-sum_b p_bi log(p_bi)}.
If \code{method="Margalef"}, the diversity index is that of Margalef (1972): \eqn{(\sum_b L_b-1)/log(\sum_b L_b a_{bi})}{(\sum_b L_b-1)/log(sum_b L_b a_bi)}.
If \code{method="Menhinick"}, the diversity index is that of Menhinick (1964): \eqn{\sum_b L_b/\sqrt{sum_b L_b a_bi}}{sum_b L_b/sqrt(sum_b L_b a_bi)}.
If \code{method="McIntosh"}, the diversity index is that of McIntosh (1967): \eqn{(\sum_b L_b a_{bi}-\sqrt{\sum_b L_b a_{bi}^2})/(\sum_b L_b a_{bi}-\sqrt{\sum_b L_b a_{bi}})}{(sum_b L_b a_bi-sqrt(sum_b L_b a_bi^2))/(sum_b L_b a_bi-sqrt(sum_b L_b a_bi))}.
If one of the strings is "full", then all indices are calculated.

\bold{For the indices relying on relative abundances to be valid, each species must support at least one feature. If this is not the case and one of these indices has to be calculated, the phylogenetic tree is re-scaled so that the shortest distance from tip to root is equal to 1. This means that the scale of the branch lengths will be changed, which will impact the calculation of feature richness (method="richness") and the indices of Margalef, Menhinick and McIntosh.}

In the rare cases where this scaling will be necessary, a better option is that the scaling be done by the user itself (prior to the use of the evodiv function). Indeed, this will unable them to choose an appropriate scaling. For example, if branch lengths are expressed as billion years of evolution and the sum of branch lengths from tip to root is lower than one billion. Then expressing branch lengths in million years of evolution may solve the problem, leading to sum of branch lengths from tip to root higher than 1 million years.  

If the phylogenetic tree has no branch lengths, all branches are set to a length of 1.
}

\references{
Gini, C. (1912) Variabilita e mutabilita. Studi economicoaguridici delle facoltta di giurizprudenza dell, Universite di Cagliari III, Parte II.

Magurran, A.E. (2004) Measuring biological diversity. Blackwell Publishing, Oxford, U.K.

Margalef, R. (1972) Homage to Evelyn Hutchinson, or why is there an upper limit to diversity? \emph{Transactions of the Connecticut Academy of Arts and Sciences}, \bold{44}, 211--235.

McIntosh, R.P. (1967) An index of diversity and the relation of certain conepts to diversity. \emph{Ecology}, \bold{48}, 392--404. 

Menhinick, E.F. (1964) A Comparison of Some Species-Individuals Diversity Indices Applied to Samples ofField Insects. \emph{Ecology}, \bold{45}, 859--861.

Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 

Pavoine, S., Ricotta, C. (2019) A simple translation from indices of species diversity to indices of phylogenetic diversity. \emph{Ecological Indicators}, \bold{101}, 552--561.

Shannon, C.E. (1948) A mathematical theory of communication. \emph{Bell System technical journal}, \bold{27}, 379--423, 623--656.

Simpson, E.H. (1949) Measurement of diversity. \emph{Nature}, \bold{163}, 688.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evodivparam}}, \code{\link{speciesdiv}}
}
\examples{
\dontrun{
if(require(ape)){
data(batcomm)
phy <- read.tree(text=batcomm$tre)
ab <- batcomm$ab[,phy$tip.label]
evodiv(phy, ab)
}
}
}
\keyword{models}
