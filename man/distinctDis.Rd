\name{distinctDis}
\alias{distinctDis}
\title{
Dissimilarity-based Species' Originality 
}
\description{
The function calculates five indices of species' originality.
}
\usage{
distinctDis(dis, method = c("Rb", "AV", "FV", "NN", "Dstar", "full"), 
    palpha = 0, standardized = FALSE)
}
\arguments{
  \item{dis}{
an object of class \code{dist} containing pair-wise (functional or phylogenetic) dissimilarities between species.
}
  \item{method}{
a vector of strings. Possible values are \code{"Rb"}, \code{"AV"}, 
\code{"FV"}, \code{"NN"}, \code{"Dstar"} and \code{"full"}. 
\code{"Rb"} is for Pavoine et al. (2017) index Rb; \code{"AV"} is for AV, the 
average dissimilarity between a species and all others in a set (Eiswerth and 
Haney 1992; Ricotta 2004); \code{"FV"} is for FV, the average dissimilarity 
between a species and any other (including the focal species itself) (Schmera 
et al. 2009); \code{"NN"} is for the minimum dissimilarity between a species 
and any other (the dissimilarity to its Nearest Neighbor) (Pavoine et al. 
2017);  \code{Dstar} is the parametric indices \eqn{^{\alpha}D*} developed by Pavoine and Ricotta (2021, Appendix B); \code{"full"} returns all indices.} 
  \item{palpha}{a numeric which provides the value of parameter \eqn{\alpha} in index \eqn{^{\alpha}D*} developed by Pavoine and Ricotta (2021, Appendix B).}
  \item{standardized}{
a logical. If \code{TRUE}, the vector of originalities is divided by its sum (transforming absolute originalities into relative originalities).
}
}
\value{
A data frame with species as rows and originality indices as columns.
}
\references{
Eiswerth, M.E. and Haney, J.C. (1992) Allocating conservation expenditures: accounting for inter-species genetic distinctiveness. \emph{Ecological Economics}, \bold{5}, 235--249.

Pavoine, S., Bonsall, M.B., Dupaix, A., Jacob, U., Ricotta, C. (2017) From phylogenetic to functional originality: guide through indices and new developments. \emph{Ecological Indicators}, \bold{82}, 196--205.

Pavoine, S., Ricotta, C. (2021) On the relationships between rarity, uniqueness, distinctiveness, originality and functional/phylogenetic diversity. \emph{BiorXiv}. \doi{10.1101/2021.08.09.455640}

Ricotta, C. (2004) A parametric diversity measure combining the relative abundances and taxonomic distinctiveness of species. \emph{Diversity and Distributions}, \bold{10}, 143--146.

Schmera, D., Podani, J., Eros, T. (2009) Measuring the contribution of community members to functional diversity. \emph{Oikos}, \bold{118}, 961--971.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{distinctTopo}}, \code{\link{distinctTree}}, \code{\link{distinctUltra}}
}
\examples{
e <- rlnorm(10)
e <- sort(e)
names(e) <- paste("s", 1:10, sep="")
d <- dist(e)
barplot(e)

D <- distinctDis(d, standardized = TRUE)

par(mfrow=c(4,2))
plot(e, D[,1], xlab="trait", ylab="Rb")
plot(e, D[,2], xlab="trait", ylab="AV")
plot(e, D[,3], xlab="trait", ylab="FV")
plot(e, D[,4], xlab="trait", ylab="NN")

plot(D[,1], D[,2], xlab="Rb", ylab="AV")
plot(D[,1], D[,3], xlab="Rb", ylab="FV")
plot(D[,2], D[,3], xlab="AV", ylab="FV")
plot(D[,2], D[,4], xlab="AV", ylab="NN")
par(mfrow=c(1,1))
}
\keyword{models}

