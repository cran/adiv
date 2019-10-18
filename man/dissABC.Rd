\name{dissABC}
\alias{dissABC}
\title{
Phylogenetic and Functional Similarity between Communities
}
\description{
Coefficients of similarity between communities that rely on the presence/absence of species are generally based on various combinations of the matching/mismatching components of the classical 2 x 2 contingency table. Three of these components are: a=the number of species shared by the two communities; b=the number of species in the first community that are not in the second; c=the number of species in the second community that are not in the first. These coefficients are extended in \code{dissABC} to include phylogenetic or functional information on species (Ricotta and Pavoine 2015).}
\usage{
dissABC(comm, dis, option = 1:4, method = c("J", "S", "O", "K", "SS","Si"))
}
\arguments{
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and relative abundance or absolute abundance as entries. Column labels (species names) should be assigned as in the object \code{dis}.
}
  \item{dis}{a matrix (or data frame) of (phylogenetic or functional) dissimilarities among species rescaled in the range [0, 1] or an object of class \code{dist} containing these dissimilarities [obtained by functions like \code{\link{dsimFun}} in this package adiv, \code{vegdist} in package vegan, \code{gowdis} in package FD, or \code{dist.ktab} in package ade4 for functional dissimilarities, or functions like \code{\link{dsimTree}} in this package adiv, \code{cophenetic.phylo} in package ape or \code{distTips} in package adephylo for phylogenetic dissimilarities]. See also function \code{\link{dsimTaxo}} for taxonomic data. Note that dissimilarities among species need first to be rescaled in the range [0, 1]. If the dissimilarities are outside the range 0-1, a warning message is displayed and each dissimilarity is divided by the maximum over all pair-wise dissimilarities.}
  \item{option}{
a numeric, either 1, 2, 3, or 4 (if several values are given only the first one is considered). See details.
}
  \item{method}{
a character or string, either \code{"J"}, \code{"S"}, \code{"O"}, \code{"K"}, \code{"SS"}, or \code{"Si"} (if several values are given only the first one is considered). See details.
}
}
\details{
To obtain the dissimilarities among plots, one needs to choose the equations to be used for the (phylogenetic or functional) components A, B, and C thanks to argument \code{option} and the way the components will be combined, thanks to argument \code{method}.

Let \eqn{\mathbf{D}=(d_{ij})}{D=(d_ij)} a matrix of (functional, morphological or phylogenetic) dissimilarities between pairs of species with \eqn{d_{ij} = d_{ji}}{d_ij = d_ji} and \eqn{d_{ii} = 0}{d_ii = 0}. If the dissimilarity coefficient \emph{d} is in the range [0, 1], it is possible to define a corresponding similarity coefficient \emph{s} as the complement of \emph{d}: \emph{s} = 1 - \emph{d}. Let \eqn{x_{ik}}{x_ik} the abundance of species \emph{i} in community \emph{k}. \emph{S(kh)} is the number of species in the pooled communities \emph{k} and \emph{h} (i.e. the species for which \eqn{min\{x_{ik}, x_{ih}\} > 0}{min(x_ik, x_ih) > 0}). The (absolute) abundance of species similar to \emph{i} in plot \emph{k} is
\deqn{Z_{ik}=\sum_{j=1}^{S(kh)}x_{jk}s_{ij}}{Z_ik=sum_j^S(kh) x_jk s_ij}. 

If \code{option=1}, equations 6-8 of the main text of Ricotta and Pavoine (2015) are used for calculating components A, B, C:
\deqn{A=\sum_{i=1}^{S(kh)}min\{Z_{ik}, Z_{ih}\}}{A=sum_i^S(kh) min(Z_ik, Z_ih)}
\deqn{B=\sum_{i=1}^{S(kh)}(max\{Z_{ik}, Z_{ih}\}-Z_{ih})}{B=sum_i^S(kh) [max(Z_ik, Z_ih)-Z_ih]}
\deqn{C=\sum_{i=1}^{S(kh)}(max\{Z_{ik}, Z_{ih}\}-Z_{ik})}{C=sum_i^S(kh) [max(Z_ik, Z_ih)-Z_ik]}

If \code{option=2}, equations A1-A3 from Appendix S1 of Ricotta and Pavoine (2015) are used. 

If \code{option=3}, equations A5-A7 from Appendix S1 of Ricotta and Pavoine (2015) are used. 

If \code{option=4}, equations A10-A12 from Appendix S1 of Ricotta and Pavoine (2015) are used.

If \code{method="J"}=the Jaccard index is used: 
\deqn{\frac{A}{A+B+C}}{A/(A+B+C)}

If \code{method="S"}=the Sorensen index is used: 
\deqn{\frac{2A}{2A+B+C}}{2A/(2A+B+C)}

If \code{method="O"}=the Ochiai index is used: \deqn{A/(\sqrt{A+B}\sqrt{A+C})}{A/[sqrt(A+B)sqrt(A+C)]}

If \code{method="K"}, the Kulczynski index is used:
\deqn{\frac{1}{2}\left(\frac{A}{A+B}+\frac{A}{A+C}\right)}{[A/(A+B)+A/(A+C)]/2}

If \code{method="SS"}, the Sokal-Sneath index is used:
\deqn{\frac{A}{A+2B+2C}}{A/(A+2B+2C)}

If \code{method="Si"}, the Simpson index is used:
\deqn{\frac{A}{A+min(B,C)}}{A/(A+min(B,2))}
}
\value{
Function \code{dissABC} returns a matrix with the values of the proposed similarities among communities based on interspecies resemblances.
}
\references{
Ricotta, C. and Pavoine, S. (2015) Measuring similarity among plots including similarity among species: an extension of traditional approaches. \emph{Journal of Vegetation Science}, \bold{26}, 1061--1067.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{discomQE}}, \code{\link{evodiss}}
}
\examples{
data(RP15JVS)
dissABC(RP15JVS$ab, RP15JVS$D1, method="J", option=1)
J <- as.matrix(dissABC(RP15JVS$ab, RP15JVS$D1, method="J", option=1))[, 1]
SS <- as.matrix(dissABC(RP15JVS$ab, RP15JVS$D1, method="SS", option=1))[, 1]
S <- as.matrix(dissABC(RP15JVS$ab, RP15JVS$D1, method="S", option=1))[, 1]
O <- as.matrix(dissABC(RP15JVS$ab, RP15JVS$D1, method="O", option=1))[, 1]
K <- as.matrix(dissABC(RP15JVS$ab, RP15JVS$D1, method="K", option=1))[, 1]
plot(1:9, J, 
xlab="Number of the plots which plot 1 is compared to",
 ylab="Similarity", type="b", ylim=c(0,1), pch=18)
lines(1:9, SS, type="b", pch=15)
lines(1:9, S, type="b", pch=17)
lines(1:9, O, type="b", pch=12)
lines(1:9, K, type="b", pch=1)
legend("bottomleft", 
c("Jaccard","Sokal-Sneath","Sorensen","Ochiai","Kulczynski"), 
pch=c(18,15,17,12,1), lty=1)
}
\keyword{models}

