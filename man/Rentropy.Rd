\name{Rentropy}
\alias{Rentropy}
\title{
Pavoine et al. (2017) functional or phylogenetic entropy
}
\description{
Function \code{Rentropy} calculates Pavoine et al. (2017) functional or phylogenetic R entropy within communities; this index is closely related to Rao's quadratic entropy 
}
\usage{
Rentropy(comm, dis = NULL, scale = FALSE)
}
\arguments{
  \item{comm}{a data frame or a matrix with communities as rows and species as columns. Entries are abundances of species within sites.}
  \item{dis}{
either \code{NULL} or an object of class \code{dist} that contains the (functional or phylogenetic) distances among species. If \code{NULL} species are said equidistant with a distance of 1 between any two species and a distance of 0 between a species and itself. 
}
  \item{scale}{
a logical value indicating whether or not the diversity coefficient should be scaled by its maximal value over all species abundance distributions.
}
}
\details{
The definition of the R entropy is:
\deqn{R(\mathbf{p}_i,\mathbf{D})=\sum_{k=1}^S\sum_{l=1}^S \sqrt{p_{k|i}}\sqrt{p_{k|j}}d_{kl}}{R(p_i,D)=sum_k,l sqrt(p_k|i) sqrt(p_k|j) d_kl}
where \eqn{\mathbf{p}_i=(p_{1|i}, ..., p_{k|i}, ..., p_{S|i})}{p_i=(p_1|i, ..., p_k|i, ..., p_S|i)} is the vector of relative species abundance within community \emph{i}; \emph{S} is  the  number  of  species; \eqn{\mathbf{D}=(d_{kl})}{D=(d_kl)} is the  matrix  of  (phylogenetic  or functional)  dissimilarities  among  species,  and \eqn{d_{kl}}{d_kl} is the (phylogenetic or functional) dissimilarity between species 
\emph{k} and \emph{l}.  
}
\value{
Function \code{Rentropy} returns a data frame with communities as rows and the R entropy within communities as columns.
}
\references{
Pavoine, S., Bonsall, M.B., Dupaix, A., Jacob, U., Ricotta, C. (2017) From phylogenetic to functional originality: guide through indices and new developments. \emph{Ecological Indicators}, \bold{82}, 196--205.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
\dontrun{
if(require(ade4)){

data(ecomor, package="ade4")
dtaxo <- dist.taxo(ecomor$taxo)
bird.R <- Rentropy(t(ecomor$habitat), dtaxo^2/2)
dotchart(bird.R$diversity, labels = rownames(bird.R))

}
}

\dontrun{
if(require(ape)){

data(batcomm)
phy <- read.tree(text=batcomm$tre)
dphy <- as.dist(cophenetic(phy))/2
ab <- batcomm$ab[, phy$tip.label]
bat.R <- Rentropy(ab, dphy)
dotchart(bat.R$diversity, labels = rownames(bat.R), xlab = "Entropy")

}
}
}
\keyword{models}
