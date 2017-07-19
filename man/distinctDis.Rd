\name{distinctDis}
\alias{distinctDis}
\title{
Dissimilarity-based Species' Originality 
}
\description{
The function calculates three indices of species' originality.
}
\usage{
distinctDis(dis, method = 1:3, standardized = FALSE)
}
\arguments{
  \item{dis}{
an object of class \code{dist} containing pair-wise (functional or phylogenetic) dissimilarities between species.
}
  \item{method}{
a numeric or a vector of numerics. Possible values are 1, 2 and 3. \code{1} is for Pavoine et al. (2017) index Rb; \code{2} is for AV, the average dissimilarity between a species and all others in a set (Eiswerth and Haney 1992; Ricotta 2004); \code{3} is for FV, the average dissimilarity between a species and any other (including the focal species itself) (Schmera et al. 2009).
}
  \item{standardized}{
a logical. If \code{TRUE}, the vector of originalities is divided by its sum (transforming absolute originalities into relative originalities).
}
}
\value{
A data frame with species as rows and originality indices as columns.
}
\references{
Eiswerth, M.E. and Haney, J.C. (1992) Allocating conservation expenditures: accounting for inter-species genetic distinctiveness. \emph{Ecological Economics}, \bold{5}, 235--249.

Ricotta, C. (2004) A parametric diversity measure combining the relative abundances and taxonomic distinctiveness of species. \emph{Diversity and Distributions}, \bold{10}, 143--146.

Schmera, D., Podani, J., Eros, T. (2009) Measuring the contribution of community members to functional diversity. \emph{Oikos}, \bold{118}, 961--971.

Pavoine, S., Bonsall, M.B., Dupaix, A., Jacob, U., Ricotta, C. (2017) From phylogenetic to functional originality: guide through indices and new developments. \code{Ecological Indicators}, \bold{82}, 196--205.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
See Also \code{\link{distinctTopo}}, \code{\link{distinctTree}}, \code{\link{distinctUltra}}
}
\examples{
e <- rlnorm(10)
e <- sort(e)
names(e) <- paste("s", 1:10, sep="")
d <- dist(e)
barplot(e)

D <- distinctDis(d, standardized = TRUE)

par(mfrow=c(3,1))
plot(e, D[,1], xlab="trait", ylab="Rb")
plot(e, D[,2], xlab="trait", ylab="AV")
plot(e, D[,3], xlab="trait", ylab="FV")

plot(D[,1], D[,2], xlab="Rb", ylab="AV")
plot(D[,1], D[,3], xlab="Rb", ylab="FV")
plot(D[,2], D[,3], xlab="AV", ylab="FV")
par(mfrow=c(1,1))
}
\keyword{models}

