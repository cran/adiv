\name{QE}
\alias{QE}
\alias{discomQE}
\title{
Quadratic Entropy
}
\description{
Function \code{QE} calculates Rao's quadratic entropy within communities

Function \code{discomQE} calculates Rao's dissimilarities between communities
}
\usage{
QE(comm, dis = NULL, formula = c("QE", "EDI"), scale = FALSE)

discomQE(comm, dis = NULL, structures = NULL, formula = c("QE", "EDI"))
}
\arguments{
  \item{comm}{a data frame or a matrix with communities as rows and species as columns. Entries are abundances of species within sites}
  \item{dis}{
either \code{NULL} or an object of class \code{dist} that contains the (functional or phylogenetic) distances among species. If \code{NULL} the Gini-Simpson index is used. 
}
  \item{formula}{
either \code{"QE"} (default) or \code{"EDI"}. See details.
}
  \item{scale}{
a logical value indicating whether or not the diversity coefficient should be scaled by its maximal value over all species abundance distributions.
}
  \item{structures}{
either NULL or a data frame that contains, in the \emph{j}th row and the \emph{k}th column, the name of the group of level \emph{k} to which the \emph{j}th community belongs. Communities in \code{structures} should be in the same order as in \code{comm}. 
}
}
\details{
If \code{formula = "QE"}, the definition of the quadratic entropy is:
\deqn{QE(\mathbf{p}_i,\mathbf{D})=\sum_{k=1}^S\sum_{l=1}^S p_{k|i}p_{k|j}d_{kl}}{QE(p_i,D)=sum_k,l p_k|i p_k|j d_kl}
where \eqn{\mathbf{p}_i=(p_{1|i}, ..., p_{k|i}, ..., p_{S|i})}{p_i=(p_1|i, ..., p_k|i, ..., p_S|i)} is the vector of relative species abundance within community \emph{i}; \emph{S} is  the  number  of  species; \eqn{\mathbf{D}=(d_{kl})}{D=(d_kl)} is the  matrix  of  (phylogenetic  or functional)  dissimilarities  among  species,  and \eqn{d_{kl}}{d_kl} is the (phylogenetic or functional) dissimilarity between species 
\emph{k} and \emph{l}. For the calculations of dissimilarities between communities see the description of the apportionment of quadratic entropy in Pavoine et al. (2016) and references therein. 

If \code{formula = "EDI"}, the definition of the quadratic entropy is: 
\deqn{EDI(\mathbf{p}_i,\mathbf{D})=\sum_{k=1}^S\sum_{l=1}^S p_{k|i}p_{k|j}\frac{d_{kl}^2}{2}}{EDI(p_i,D)=sum_k,l p_k|i p_k|j (d_kl^2)/2}
EDI stands for the Euclidean Diversity Index of Champely and Chessel (2002) (equation 3 in Pavoine et al. 2004). If EDI is used, the dissimilarities calculated by \code{discomQE} are obtained from the same scheme (equation 4 in Pavoine et al. 2004).
}
\value{
Function \code{QE} returns a data frame with communities as rows and the diversity within communities as columns.

If \code{structures} is \code{NULL}, function \code{discomQE} returns an object of class \code{dist}. Otherwise it returns a list of objects of class \code{dist}.

}
\references{
Gini, C. (1912) \emph{Variabilita e mutabilita}. Universite di Cagliari III, Parte II.

Simpson, E.H. (1949) Measurement of diversity. \emph{Nature}, \bold{163}, 688. 

Rao, C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. \emph{Theoretical Population Biology}, \bold{21}, 24--43.

Champely, S. and Chessel, D. (2002) Measuring biological diversity using Euclideanmetrics. \emph{Environmental and Ecological Statistics}, \bold{9}, 167--177.

Pavoine, S., Dufour, A.B., Chessel, D. (2004) From dissimilarities among species to dissimilarities among communities: a double principal coordinate analysis. \emph{Journal of Theoretical Biology}, \bold{228:}, 523--537.

Pavoine, S., Marcon, E., Ricotta, C. (2016) "Equivalent numbers" for species, phylogenetic, or functional diversity in a nested hierarchy of multiple scales. \emph{Methods in Ecology and Evolution}, \bold{7}, 1152--1163.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
data(ecomor, package="ade4")
dtaxo <- dist.taxo(ecomor$taxo)
QE(t(ecomor$habitat), dtaxo, formula="EDI")
QE(t(ecomor$habitat), dtaxo^2/2, formula="QE")
table.value(as.matrix(discomQE(t(ecomor$habitat), dtaxo, formula="EDI")))
EDIcom <- discomQE(t(ecomor$habitat), dtaxo, formula="EDI")
QEcom <- discomQE(t(ecomor$habitat), dtaxo^2/2, formula="QE")
QEcom 
EDIcom^2/2

bird.QE <- QE(t(ecomor$habitat), dtaxo, formula="EDI")
dotchart(bird.QE$diversity, labels = rownames(bird.QE))

data(humDNAm, package="ade4")
QE(t(humDNAm$samples), humDNAm$distances/2, formula="QE")
QE(t(humDNAm$samples), sqrt(humDNAm$distances), formula="EDI")
QEhumDNA.dist <- discomQE(t(humDNAm$samples), 
humDNAm$distances/2, humDNAm$structures)
is.euclid(QEhumDNA.dist$communities)
is.euclid(QEhumDNA.dist$regions)

EDIhumDNA.dist <- discomQE(t(humDNAm$samples), 
sqrt(humDNAm$distances), humDNAm$structures, formula="EDI")
is.euclid(EDIhumDNA.dist$communities)
is.euclid(EDIhumDNA.dist$regions)

QEhumDNA.dist$communities
EDIhumDNA.dist$communities^2/2

hum.QE <- QE(t(humDNAm$samples), humDNAm$distances/2, formula="QE")
dotchart(hum.QE$diversity, labels = rownames(hum.QE))
}
\keyword{models}
