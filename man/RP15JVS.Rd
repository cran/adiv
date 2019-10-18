\name{RP15JVS}
\alias{RP15JVS}
\docType{data}
\title{
Theoretical Data Set used in Ricotta and Pavoine (2015) in Journal of Vegetation Science
}
\description{
The data set was used in Ricotta and Pavoine (2015) to illustrate the relationships between the coefficients of similarity between communities that they developed. These data represent an artificial ecological gradient. 
}
\usage{data("RP15JVS")}
\format{
The format is a list of four objects:

\code{ab}, a data frame with communities as rows, species as columns, and number of indivuduals as entries.

\code{D1}, a data frame of pairwise dissimilarities between species. In \code{D1}, dissimilarities reflect the species ecological differences along the artificial gradient of table \code{ab}: interspecies dissimilarities were set roughly proportional to the distance between the locations of the species optima and to the difference between their optimal abundances (see Ricotta and Pavoine 2015 for details).

\code{D2}, a data frame of pairwise dissimilarities between species. In \code{D2}, dissimilarities were randomly assigned using an even distribution.

\code{D3}, a data frame of pairwise dissimilarities between species. In \code{D3}, dissimilarities were uniformly set to 2/3.
}
\source{
Appendixes of Ricotta and Pavoine (2015)
}
\references{
Ricotta, C. and Pavoine, S. (2015) Measuring similarity among plots including similarity among species: An extension of traditional approaches. \emph{Journal of Vegetation Science}, \bold{26}, 1061--1067
}
\examples{
data(RP15JVS)
dissABC(RP15JVS$ab, RP15JVS$D1, method="J", option=1)
}
\keyword{datasets}
