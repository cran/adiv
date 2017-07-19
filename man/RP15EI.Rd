\name{RP15EI}
\alias{RP15EI}
\docType{data}
\title{
Theoretical Data Set used in Ricotta and Pavoine (2015) in Ecological Indicators
}
\description{
The data set was used in Ricotta and Pavoine (2015) to show how the coefficients of multi-site dissimilarity that they developed can be applied to data.  
}
\usage{data("RP15EI")}
\format{
The format is a list of 4 data frames. Each data frame gives the presence (1) and absence (0) of 10 species (columns, S1-S10) in 8 plots (rows P1-P8). These data are characterized by different levels of nestedness and turnover in species compositions between sites.

\code{M1}: Intermediate nestedness and turnover.

\code{M2}: Perfectly nested pattern.

\code{M3}: Very high species turnover.

\code{M4}: Random configuration.
}
\source{
Figure 1 of Ricotta and Pavoine (2015)
}
\references{
Ricotta, C. and Pavoine, S. (2015) A multiple-site dissimilarity measure for species presence/absence data and its relationship with nestedness and turnover. \emph{Ecological Indicators}, \bold{54},203--206.
}
\examples{
data(RP15EI)
betastatjac(RP15EI$M1)
}
\keyword{datasets}
