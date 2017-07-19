\name{RDMCCP16}
\alias{RDMCCP16}
\docType{data}
\title{
Theoretical Data Set used in Ricotta et al. (2016)
}
\description{
Hypothetical communities C1-C9 composed of nine species
S1-S9 with varying abundances divided into three groups of three species (S1-S3, S4-S6 and S7-S9) (say, legumes, herbs and forbs). For simplicity, all species within the same group are functionally identical to each other, while two species belonging to different groups are always maximally dissimilar (Ricotta et al. 2016).}
\usage{data("RDMCCP16")}
\format{
A list of two objects:

\code{ab}, a matrix with communities as rows, species as columns and abundance values as entries.

\code{dis}, an object of class \code{dist} that contains theoretical dissimilarities between species.
}
\source{
Table 1 in Ricotta et al. (2016)
}
\references{
Ricotta, C., de Bello, F., Moretti, M., Caccianiga, M., Cerabolini, B.E., Pavoine, S. (2016). Measuring the functional redundancy of biological communities: A quantitative guide. \emph{Methods in Ecology and Evolution}, \bold{7}, 1386--1395.}
\examples{
data(RDMCCP16)
uniqueness(RDMCCP16$ab, as.dist(RDMCCP16$dis))
}
\keyword{datasets}
