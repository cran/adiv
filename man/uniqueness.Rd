\name{uniqueness}
\alias{uniqueness}
\title{
Functional Uniqueness and Functional Redundancy of Biological Communities
}
\description{
Hypothetical communities C1-C9 composed of nine species
S1-S9 with varying abundances divided into three groups of three species each displayed with different shades (S1-S3, S4-S6 and S7-S9). All species within the same group are functionally identical to each other, while two species belonging to different groups are always maximally dissimilar.}
\usage{
uniqueness(comm, dis, tol = 1e-08, abundance = TRUE)
}
\arguments{
  \item{comm}{
a matrix or a data frame containing the abundance or incidence (0/1) of species in communities (or plots). Columns are species and communities are rows.
}
  \item{dis}{
an object of class \code{dist} containing the functional distances among species}
  \item{tol}{
a tolerance threshold (a value between -\code{tol} and \code{tol} is considered as null)
}
  \item{abundance}{
a logical. If \code{TRUE}, abundance data are used when available; if \code{FALSE}, incidence (0/1) data are used.}
}
\value{
The function uniqueness returns a list of three data frames:

\code{kbar}: this first data frame gives values for Ricotta et al. coefficient \eqn{\bar{K}_i}{K_i}'s per species (rows) and community (columns).

\code{V}: this second data frame gives values for Ricotta et al. coefficient \eqn{V_i}'s per species (rows) and community (columns).

\code{red}: this third data set gives values for Ricotta et al. coefficients \code{N} (species richness), \code{Q} (quadratic diversity), \code{D} (Simpson diversity), and \code{U} (uniqueness) per community; in this third data frame, coefficients are columns and communities are rows, the coefficients are thus calculated per community only.
}
\references{
Ricotta, C, de Bello, F, Moretti, M, Caccianiga, M, Cerabolini, BE, Pavoine, S (2016). Measuring the functional redundancy of biological communities: A quantitative guide. \emph{Methods in Ecology and Evolution}, \bold{7}, 1386--1395.}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
See Also \code{\link{QE}}
}
\examples{
data(RDMCCP16)
uniqueness(RDMCCP16$ab, as.dist(RDMCCP16$dis))
}
\keyword{models}