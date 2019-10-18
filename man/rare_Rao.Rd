\name{rare_Rao}
\alias{rare_Rao}
\title{
Functional Rarefaction for Species Abundance Data
}
\description{
The function \code{Rare_Rao} performs distance-based rarefaction curves using species abundance data. It finds the expected functional diversity (if functional distances between species are used) as a function of the sampling effort. Two approaches are available: an analytical solution, a resampling approach.}
\usage{
rare_Rao(comm, dis, sim = TRUE, resampling = 999, formula = c("QE", "EDI"))
}
\arguments{
  \item{comm}{
a data frame or a matrix with samples as rows, species as columns, and abundance or frequency as entries. If presences/absences (1/0) are given, the relative abundance of a given species in a community of S species will be considered equal to 1/S.
}
  \item{dis}{
an object of class \code{dist} containing pairwise distances among species. The distance matrix should be squared Euclidean or simply Euclidean (see definition and details in Appendix S1 of Ricotta et al. (2012)). Alternatively, if \code{dis} is \code{NULL}, the Gini-Simpson index is used.
}
  \item{sim}{
a logical; if \code{TRUE}, the resampling approach is used; if \code{FALSE}, the analytical solution is given. 
}
  \item{resampling}{
a numeric; number of times data are resampled to calculate the mean functional rarefaction curve (used if \code{sim=TRUE}).
}
  \item{formula}{
either \code{"QE"} (default) or \code{"EDI"}. See details.
}
}
\details{
If \code{formula = "QE"}, the definition of the quadratic entropy is:
\deqn{QE(\mathbf{p}_i,\mathbf{D})=\sum_{k=1}^S\sum_{l=1}^S p_{k|i}p_{k|j}d_{kl}}{QE(p_i,D)=sum_k,l p_k|i p_k|j d_kl}
where \eqn{\mathbf{p}_i=(p_{1|i}, ..., p_{k|i}, ..., p_{S|i})}{p_i=(p_1|i, ..., p_k|i, ..., p_S|i)} is the vector of relative species abundance within sample \emph{i}; \emph{S} is  the  number  of  species; \eqn{\mathbf{D}=(d_{kl})}{D=(d_kl)} is the  matrix  of  (phylogenetic  or functional)  dissimilarities  among  species,  and \eqn{d_{kl}}{d_kl} is the (phylogenetic or functional) dissimilarity between species 
\emph{k} and \emph{l}.  

If \code{formula = "EDI"}, the definition of the quadratic entropy is: 
\deqn{EDI(\mathbf{p}_i,\mathbf{D})=\sum_{k=1}^S\sum_{l=1}^S p_{k|i}p_{k|j}\frac{d_{kl}^2}{2}}{EDI(p_i,D)=sum_k,l p_k|i p_k|j (d_kl^2)/2}
EDI stands for the Euclidean Diversity Index of Champely and Chessel (2002) (equation 3 in Pavoine et al. 2004). 

In both cases, if \code{dis = NULL}, the quadratic entropy is equal to Gini-Simpson entropy:
\deqn{H_{GS}(\mathbf{p}_i)=1 - \sum_{k=1}^S p_{k|i}^2}{H_GS(p_i)=1 - sum_k (p_k|i)^2}
}
\value{
If \code{sim = TRUE}, the function returns a data frame containing the Expected Rao Quadratic Entropy (column 'ExpRao'), the limits of the 95\% Confidence Interval (columns 'LeftIC' and 'RightIC') for each subsample dimension (M) out of the total set of samples (N). If \code{sim = FALSE}, the function returns a data frame containing the analytical solution for the Expected Rao Quadratic Entropy (column 'ExpRao') for each subsample dimension (M) out of the total set of samples (N).
}
\references{
Ricotta, C., Pavoine, S., Bacaro, G., Acosta, A. (2012) Functional rarefaction for species abundance data. \emph{Methods in Ecology and Evolution}, \bold{3}, 519--525.

Champely, S. and Chessel, D. (2002) Measuring biological diversity using Euclideanmetrics. \emph{Environmental and Ecological Statistics}, \bold{9}, 167--177.

Pavoine, S., Dufour, A.B., Chessel, D. (2004) From dissimilarities among species to dissimilarities among communities: a double principal coordinate analysis. \emph{Journal of Theoretical Biology}, \bold{228}, 523--537.
}
\author{
Giovanni Bacaro
and
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{QE}}
}
\examples{
\dontrun{
if(require(ade4)){

data(aviurba, package="ade4")

# Trait-based distances between bird species:
distances<-dist.ktab(ktab.list.df(list(aviurba$traits)), type = "N")
# The distances should be squared Euclidean; 
# note that Euclidean distances can be used 
# as they also are squared Euclidean.

# Species abundances in sites
abundances<- aviurba$fau

# Rarefaction of functional diversity
rare_Rao(abundances, distances, sim = TRUE, resampling = 100)
rare_Rao(abundances, distances, sim = FALSE)

}
}
}
\keyword{models}
