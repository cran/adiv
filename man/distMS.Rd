\name{distMS}
\alias{distMS}
\title{
Marczewski-Steinhaus Coefficient
}
\description{
The function calculates the Marczewski-Steinhaus coefficient of dissimilarity between pairs of entities (e.g. communities)}
\usage{
distMS(comm, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{comm}{
a data frame or a matrix of nonnegative values (e.g. abundance of species (columns) within communities (rows) to obtain dissimilarities between communities).
}
  \item{diag}{
a logical value indicating whether the diagonal of the distance matrix should be printed by function \code{print.dist}.
}
  \item{upper}{
a logical value indicating whether the upper triangle of the distance matrix should be printed by function \code{print.dist}.
}
}
\value{
an object of class \code{dist}
}
\references{
Orloci, L. (1978) \emph{Multivariate Analysis in Vegetation Research}. The Hague: Junk.

Legendre, P. and Legendre, L. (1998) Numerical Ecology. Amsterdam: Elsevier.

Ricotta, C., de Bello, F., Moretti, M., Caccianiga, M., Cerabolini, B.E., Pavoine, S. (2016). Measuring the functional redundancy of biological communities: a quantitative guide. \emph{Methods in Ecology and Evolution}, \bold{7}, 1386--1395.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\note{
This function is a modification of function \code{dist.quant} from library ade4 where other dissimilarity coefficients can be found.
}
\examples{
data(birdData)
distMS(birdData$fau)
}
\keyword{models}