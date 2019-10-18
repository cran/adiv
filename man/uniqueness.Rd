\name{uniqueness}
\alias{uniqueness}
\title{
Functional Uniqueness and Functional Redundancy of Biological Communities
}
\description{
The function \code{Uniqueness} calculates community-level functional uniqueness and redundancy.
}
\usage{
uniqueness(comm, dis, tol = 1e-08, abundance = TRUE)
}
\arguments{
  \item{comm}{
a matrix or a data frame containing the abundance or incidence (0/1) of species in communities (or plots). Columns are species and communities are rows.
}
  \item{dis}{
an object of class \code{dist} containing the functional distances among species. Values in dis must be bounded between 0 and 1. If they are not bounded, the function divides all values in dis by the highest observed value in dis.}
  \item{tol}{
a tolerance threshold (a value between -\code{tol} and \code{tol} is considered as null).
}
  \item{abundance}{
a logical. If \code{TRUE}, abundance data are used when available; if \code{FALSE}, incidence (0/1) data are used.}
}
\value{
The function uniqueness returns a list of three data frames:

\code{kbar}: this first data frame gives values for Ricotta et al. (2016) coefficient \eqn{\bar{K}_i}{K_i}'s per species (rows) and community (columns).

\code{V}: this second data frame gives values for Ricotta et al. (2016) coefficient \eqn{V_i}'s per species (rows) and community (columns).

\code{red}: this third data set gives values, per community, for Ricotta et al. (2016) coefficients \code{N} (species richness), \code{Q} (quadratic diversity), \code{D} (Simpson diversity), \code{U=Q/D}  (uniqueness), \code{R=1-U}  (redundancy), and Pavoine and Ricotta (2019) \code{Ustar=(1-D)/(1-Q)}  (uniqueness) and \code{Rstar=1-Ustar} (redundancy); in this third data frame, coefficients are columns and communities are rows, the coefficients are thus calculated per community only.
}
\references{
Ricotta, C., de Bello, F., Moretti, M., Caccianiga, M., Cerabolini, B.E., Pavoine, S. (2016). Measuring the functional redundancy of biological communities: A quantitative guide. \emph{Methods in Ecology and Evolution}, \bold{7}, 1386--1395.

Pavoine, S., Ricotta, C. (2019). A simple translation from indices of species diversity to indices of phylogenetic diversity. \emph{Ecological Indicators}, \bold{101}, 552--561.
} 
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{QE}}, \code{\link{treeUniqueness}}
}
\examples{
data(RDMCCP16)
uniqueness(RDMCCP16$ab, as.dist(RDMCCP16$dis))
}
\keyword{models}