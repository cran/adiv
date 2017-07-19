\name{EqRao}
\alias{EqRS}
\alias{EqRSintra}
\alias{EqRao}
\alias{wapqe}
\alias{rtestEqRS}
\alias{rtestEqRSintra}
\alias{rtestEqRao}
\alias{rtestwapqe}
\title{
Apportionment of Diversity
}
\description{
Eight functions are available.

\code{EqRS} performs the first decomposition of diversity developed in Pavoine et al. (2016) and \code{rtestEqRS} can be used for the associated permutation tests. 

\code{EqRSintra} performs the second decomposition of diversity developed in Pavoine et al. (2016) and \code{rtestEqRSintra} can be used for the associated permutation tests. 

\code{EqRao} performs the third decomposition of diversity introduced in Pavoine et al. (2016) and \code{rtestEqRao} can be used for the associated permutation tests. 
 
\code{wapqe} performs the apportionment of quadratic entropy (Rao 1986) and \code{rtestwapqe} associated permutation tests (Pavoine et al. 2016).
}
\usage{
EqRS(comm, dis = NULL, structures = NULL, 
option = c("eq", "normed1", "normed2"), 
formula = c("QE", "EDI"), tol = 1e-08)

EqRSintra(comm, dis = NULL, structures = NULL, 
option = c("eq", "normed1", "normed2"), 
formula = c("QE", "EDI"), tol = 1e-08, 
metmean = c("harmonic", "arithmetic"))

EqRao(comm, dis = NULL, structures = NULL, 
option = c("eq", "normed1", "normed2"), 
formula = c("QE", "EDI"), wopt = c("even", "speciesab"), 
tol = 1e-08, metmean = c("harmonic", "arithmetic"))

wapqe(comm, dis = NULL, structures = NULL, 
formula = c("QE", "EDI"), wopt = c("even", "speciesab"), 
tol = 1e-08)

rtestEqRS(comm, dis = NULL, structures = NULL, 
formula = c("QE", "EDI"), 
option = c("normed1", "normed2", "eq"), level = 1, nrep = 99, 
alter = c("greater", "less", "two-sided"), tol = 1e-08)

rtestEqRSintra(comm, dis = NULL, structures = NULL, 
formula = c("QE", "EDI"), 
option = c("normed1", "normed2", "eq"), level = 1, 
nrep = 99, alter = c("greater", "less", "two-sided"), 
tol = 1e-08, metmean = c("harmonic", "arithmetic"))

rtestEqRao(comm, dis = NULL, structures = NULL, 
formula = c("QE", "EDI"), 
option = c("normed1", "normed2", "eq"), 
wopt = c("even", "speciesab"), level = 1, nrep = 99, 
alter = c("greater", "less", "two-sided"), tol = 1e-08, 
metmean = c("harmonic", "arithmetic"))

rtestwapqe(comm, dis = NULL, structures = NULL, 
formula = c("QE", "EDI"), wopt = c("even", "speciesab"), 
level = 1, nrep = 99, alter = c("greater", "less", "two-sided"),
tol = 1e-08)
}
\arguments{
  \item{comm}{
a data frame or a matrix with communities as rows and species as columns. Entries are abundances of species within sites.}
  \item{dis}{
either \code{NULL} or an object of class \code{dist} that contains the (functional or phylogenetic) distances among species. If \code{NULL}, the Gini-Simpson index is used.}
  \item{structures}{
either \code{NULL} or a data frame that contains, in the \emph{j}th row and the \emph{k}th column, the name of the group of level \emph{k} to which the \emph{j}th community belongs. Communities in \code{structures} should be in the same order as in \code{comm}. See details.}
  \item{option}{
\code{"eq"}, \code{"normed1"} or \code{"normed2"} (if a vector with several codes is given, only the first one is considered). \code{"eq"} means that the diversity components are given in terms of equivalent number (E) of species, sites, regions etc. \code{"normed1"} means that the normed components of diversity will be returned with formula (1 - 1 / E) / (1  - 1 / Emax);  \code{"normed2"} means that the normed components of diversity will be returned with formula (E - 1) / (Emax - 1).  For Ealpha and Egamma, Emax=S (the number of species in the data set).}
  \item{formula}{
either \code{"QE"} (default) or \code{"EDI"}. See details.
}
  \item{wopt}{
either \code{"even"} (default), \code{"speciesab"}, or a numeric vector. If \code{w="speciesab"}, then the communities will be weighted by their sum of species' abundances. If \code{w="even"}, the communities will be evenly weighted within the factors defined by the parameter \code{structures}.}
  \item{tol}{
a tolerance threshold (a value between -\code{tol} and \code{tol} is considered equal to zero)}
  \item{metmean}{
a character string: either \code{"arithmetic"} or \code{"harmonic"} (default).}
  \item{level}{
a number. The number is discarded if the parameter \code{structures} is set to \code{NULL}. If \code{structures} is different from \code{NULL} then \code{1} means test for differences among communities, within the levels of the first factor given in parameter \code{structures} (column 1), \code{2} means test for differences among levels of the first factor given in parameter \code{structures} (column 1) but within levels of the second factor given in parameter \code{structures} (column 2) (if available), etc. \code{level} should thus be between 1 and 1+the number of columns in \code{structures}. See details for concrete examples.}
  \item{nrep}{
the number of permutations.
}
  \item{alter}{
a character string specifying the alternative hypothesis, must be one of \code{"greater"} (default), \code{"less"} or \code{"two-sided"}.
}
}
\details{
If \code{formula = "QE"}, the definition of the quadratic entropy is:
\deqn{QE(\mathbf{p}_i,\mathbf{D})=\sum_{k=1}^S\sum_{l=1}^S p_{k|i}p_{k|j}d_{kl}}{QE(p_i,D)=sum_k,l p_k|i p_k|j d_kl}
where \eqn{\mathbf{p}_i=(p_{1|i}, ..., p_{k|i}, ..., p_{S|i})}{p_i=(p_1|i, ..., p_k|i, ..., p_S|i)} is the vector of relative species abundance within community \emph{i}; \emph{S} is  the  number  of  species; \eqn{\mathbf{D}=(d_{kl})}{D=(d_kl)} is the  matrix  of  (phylogenetic  or functional)  dissimilarities  among  species,  and \eqn{d_{kl}}{d_kl} is the (phylogenetic or functional) dissimilarity between species 
\emph{k} and \emph{l}.  

If \code{formula = "EDI"}, the definition of the quadratic entropy is: 
\deqn{EDI(\mathbf{p}_i,\mathbf{D})=\sum_{k=1}^S\sum_{l=1}^S p_{k|i}p_{k|j}\frac{d_{kl}^2}{2}}{EDI(p_i,D)=sum_k,l p_k|i p_k|j (d_kl^2)/2}
EDI stands for the Euclidean Diversity Index of Champely and Chessel (2002). For example, applying EqRao with the distances \code{dis}=mydis and \code{formula=QE} corresponds to applying it with \code{dis}=sqrt(2*mydis) and \code{formula=EDI}.

Parameter level permits to choose the component of 
beta diversity that will be tested for significance.

Examples are given below:

Scenario 1: Imagine that \code{structures} is \code{NULL}, then we only have a set of communities, species within them and a matrix of (functional or phylogenetic) dissimilarities among species. Then there is only one level of beta diversity which represents the average dissimilarities among communities. The functions \code{rtestEqRS}, \code{rtestEqRSintra}, \code{rtestEqRao}, and \code{rtestwapqe} will always, in that case, test for the significance of the dissimilarities among communities (whatever the value given to parameter \code{level}).

Scenario 2:
Imagine that \code{structures} is a data frame with sites as rows and only one column representing how sites are distributed among regions. Then, there are two levels of beta diversity: beta1 diversity represents the (functional or phylogenetic) dissimilarities among sites within regions; and beta2 diversity represents the (functional or phylogenetic) dissimilarities among regions. If \code{level = 1} then functions \code{rtestEqRS}, \code{rtestEqRSintra}, \code{rtestEqRao}, and \code{rtestwapqe} will test for the significance of the dissimilarities among sites within regions (beta1 diversity); in contrast, if \code{level = 2} functions \code{rtestEqRS}, \code{rtestEqRSintra}, \code{rtestEqRao}, and \code{rtestwapqe} will test for the significance of the dissimilarities among regions (beta2 diversity). As there is only one column in parameter \code{structures} and thus only two levels of diversity, \code{level} cannot be higher than 2.

Scenario 3:
Imagine now that parameter \code{comm} contains quadrats
as rows and species as columns and that \code{structures} is a data frame with quadrats as rows and three columns
representing how quadrats are organized in three nested factors: hamlets, towns, and counties. Then there would be four levels of beta diversity: beta1 diversity represents the (functional or phylogenetic) dissimilarities among quadrats within hamlets, towns, and counties; beta2 diversity represents the dissimilarities among hamlets within towns and counties; beta3 diversity represents the dissimilarities among towns within counties; and beta4 diversity represents the dissimilarities among counties. If \code{level = 1} then functions \code{rtestEqRS}, \code{rtestEqRSintra}, \code{rtestEqRao}, and \code{rtestwapqe} will test for the significance of the dissimilarities among quadrats within hamlets, towns, and counties (beta1 diversity); if \code{level = 2} functions \code{rtestEqRS}, \code{rtestEqRSintra}, \code{rtestEqRao}, and \code{rtestwapqe} will test for the significance of the dissimilarities among hamlets within towns and regions (beta2 diversity); if \code{level = 3} functions \code{rtestEqRS}, \code{rtestEqRSintra}, \code{rtestEqRao}, and \code{rtestwapqe} will test for the significance of the dissimilarities among towns within regions (beta3 diversity); if \code{level = 4} functions \code{rtestEqRS}, \code{rtestEqRSintra}, \code{rtestEqRao}, and \code{rtestwapqe} will test for the significance of the dissimilarities among regions (beta4 diversity). As there are only three columns in parameter \code{structures} and thus only four levels of diversity, \code{level} cannot be higher than 4.
Test for \code{level=1} is performed by permuting the 
abundances of each species across quadrats but within hamlets, towns, and counties. Test for \code{level=2} is performed by permuting the quadrats among hamlets within towns and counties. 
Test for \code{level=3} is performed by permuting the hamlets among towns within counties. Test for \code{level=4} is performed by permuting the towns among counties.

Other permutation schemes will be added in the future.
}
\value{
The functions \code{EqRS}, \code{EqRSintra}, \code{EqRao}, \code{wapqe} return a data frame with each component of the selected diversity decomposition. The functions \code{rtestEqRS}, \code{rtestEqRSintra}, \code{rtestEqRao}, and \code{rtestwapqe} return a list of class \code{randtest} or \code{krandtest} (classes of package ade4).
}
\references{
Pavoine, S., Marcon, E., Ricotta, C. (2016) "Equivalent numbers" for species, phylogenetic, or functional diversity in a nested hierarchy of multiple scales. \emph{Methods in Ecology and Evolution}, \bold{7}, 1152--1163.

Champely, S. and Chessel, D. (2002) Measuring biological diversity using Euclideanmetrics. \emph{Environmental and Ecological Statistics}, \bold{9}, 167--177.

Rao, C.R. (1986) Rao's axiomatization of diversity measures. In: \emph{Encyclopedia of Statistical Sciences}. Vol. 7 (eds S. Kotz and N.L. Johnson), pp. 614--617. New York: Wiley and Sons.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
data(macroloire, package="ade4")

# Taxonomic dissimilarities among species:
dTaxo <- dist.taxo(macroloire$taxo)^2/2
dTaxo <- dTaxo/max(dTaxo)
# Size-based dissimilarities among species
dSize <- dist.prop(macroloire$traits[ ,1:4], method = 2)
# Dissimilarities among species in terms of feeding categories
dFeed <- dist.prop(macroloire$traits[ ,5:11], method = 2)
# Dissimilarities among species in terms of 
# both size and feeding categories
dSF <- (dSize+dFeed)/2

# Table with sites as rows (stations), 
# species as columns and abundances as entries
ab <- as.data.frame(t(macroloire$fau))
# Table with sites as rows and one column only. 
# Entries indicate the geological region associated 
# with each site
stru <- macroloire$envir["Morphoregion"]

EqRao(ab, , stru, option="eq")
EqRao(ab, dTaxo, stru, formula = "QE", option="eq")
EqRao(ab, dSize, stru, formula = "QE", option="eq")
EqRao(ab, dFeed, stru, formula = "QE", option="eq")
EqRao(ab, dSF, stru, formula = "QE", option="eq")

EqRao(ab, , stru, option="normed2")
EqRao(ab, dTaxo, stru, formula = "QE", option="normed2")
EqRao(ab, dSize, stru, formula = "QE", option="normed2")
EqRao(ab, dFeed, stru, formula = "QE", option="normed2")
EqRao(ab, dSF, stru, formula = "QE", option="normed2")

\dontrun{
# Tests for dissimilarities among sites within regions:

rb1_GS <- rtestEqRao(ab, , stru, level=1, nrep=999, option="normed2")
rb1_GS
plot(rb1_GS)
rb1_Taxo <- rtestEqRao(ab, dTaxo, stru, formula = "QE", 
level=1, nrep=999, option="normed2")
rb1_Taxo
plot(rb1_Taxo)
rb1_Size <- rtestEqRao(ab, dSize, stru, formula = "QE", 
level=1, nrep=999, option="normed2")
rb1_Size
plot(rb1_Size)
rb1_Feed <- rtestEqRao(ab, dFeed, stru, formula = "QE", 
level=1, nrep=999, option="normed2")
rb1_Feed
plot(rb1_Feed)
rb1_SF <- rtestEqRao(ab, dSF, stru, formula = "QE", 
level=1, nrep=999, option="normed2")
rb1_SF
plot(rb1_SF)

# Tests for dissimilarities among regions:

r2_GS <- rtestEqRao(ab, , stru, level=2, nrep=999, option="normed2")
r2_GS
plot(r2_GS)
r2_Taxo <- rtestEqRao(ab, dTaxo, stru, formula = "QE", 
level=2, nrep=999, option="normed2")
r2_Taxo
plot(r2_Taxo)
r2_Size <- rtestEqRao(ab, dSize, stru, formula = "QE", 
level=2, nrep=999, option="normed2", w="even")
r2_Size
plot(r2_Size)
r2_Feed <- rtestEqRao(ab, dFeed, stru, formula = "QE", 
level=2, nrep=999, option="normed2", w="even")
r2_Feed
plot(r2_Feed)
r2_SF <- rtestEqRao(ab, dSF, stru, formula = "QE", 
level=2, nrep=999, option="normed2", w="even")
r2_SF
plot(r2_SF)
}
}
\keyword{models}
