\name{FunImbalance}
\alias{FunImbalance}
\title{Indices of functional imbalance}
\description{
The function \code{FunImbalance} calculates functional imbalance within species communities, as an indicator of the strength of interaction between
species abundances and their functional dissimilarities.
}
\usage{
FunImbalance(comm, dis, method = c("CorB", "SESB", "QB"), nrep = 10000, tol = 1e-16) 
}
\arguments{
\item{comm}{matrix or data frame with communities as rows, species as columns and abundance data as entries.}
\item{dis}{an object of class dist that provides the functional dissimilarities between species.}
\item{method}{a character string: one or a vector of "CorB", "SESB", and "QB". The strings indicate the name of the index of functional imbalance that must be used.}
\item{nrep}{an integer. The number of repetitions (permutations) to use in the calculation of SESB and QB.}
\item{tol}{a tolerance number: a value in [-tol, tol] is considered as zero.}
}
\value{
The function \code{betaUniqueness} returns a data frame with communities as rows and the selected indices of functional imbalance as columns.
}
\references{
Ricotta, C., Bacaro, G., Maccherini, S., Pavoine, S. (2022) Functional imbalance not functional evenness is the third component of community structure, \emph{Contact author for information}
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
data(RutorGlacier)
fundis <- dist(scale(RutorGlacier$Traits2[1:6]))
fundis <- fundis/max(fundis)
funImb <- FunImbalance(RutorGlacier$Abund, fundis, method = "CorB")
funImb
}
\keyword{models}

