\name{betastat}
\alias{betastatjac}
\alias{betastatsor}
\title{
Multiple-Site Dissimilarity Measure for Species Presence/Absence
Data
}
\description{
Functions \code{betastatjac} and \code{betastatsor} calculate multiple-site dissimilarity (beta diversity). The first one is derived from Jaccard coefficient of similarity and the second from Sorensen coefficient. These proposed dissimilarity indices are additively partitioned into species nestedness and turnover.
}
\usage{
betastatjac(comm)

betastatsor(comm)
}
\arguments{
  \item{comm}{
a data frame typically with communities as rows, species as columns and presence/absence (1/0) as entries. 
}
}
\value{
The two functions return a vector of 4 values:
\item{beta}{Ricotta and Pavoine (2015) \eqn{\beta^+} relative measure of additive beta diversity (multiple-site dissimilarity);}
\item{betaT}{Ricotta and Pavoine (2015) \eqn{\beta_T} contribution of species turnover to multiple-site dissimilarity;}
\item{betaN}{Ricotta and Pavoine (2015) \eqn{\beta_N} contribution of species nestedness to multiple-site dissimilarity;}
\item{sim}{Ricotta and Pavoine (2015) \eqn{\bar{\beta}^{\times}}{\beta^\Chi} relative measure of multiple-site similarity.}
}
\references{
Ricotta, C. and Pavoine, S. (2015) A multiple-site dissimilarity measure for species presence/absence data and its relationship with nestedness and turnover. \emph{Ecological Indicators}, \bold{54}, 203--206.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
data(RP15EI)
# Scripts used in Figure 1 of Ricotta and Pavoine (2015)
betastatjac(RP15EI$M1)
betastatjac(RP15EI$M2)
betastatjac(RP15EI$M3)
betastatjac(RP15EI$M4)

#see also
betastatsor(RP15EI$M1)
betastatsor(RP15EI$M2)
betastatsor(RP15EI$M3)
betastatsor(RP15EI$M4)
}
\keyword{models}
