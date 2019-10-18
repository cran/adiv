\name{RPP16EE}
\alias{RPP16EE}
\docType{data}
\title{
Theoretical Data Set used in Ricotta et al. (2016) in Ecology and Evolution
}
\description{
The data set was used in Ricotta et al. (2016) to show how the coefficients of plot-to-plot dissimilarity that they developed can be applied to data.  
}
\usage{data("RPP16EE")}
\format{
The format is a list of 2 objects:

\code{Com}, a data frame. Artificial data table composed of 15 species (S1-S15) (columns) and 9 plots (P1-P9) (rows).

\code{Dis}, a data frame with the artificial dissimilarities between species.
}
\source{
Appendixes 1 and 2 of Ricotta et al. (2016)
}
\references{
Ricotta, C., Podani, J., Pavoine, S. (2016) A family of functional dissimilarity measures for presence and absence data. \emph{Ecology and Evolution}, \bold{6}, 5383--5389.
}
\examples{
data(RPP16EE)
RPP16EE$Com
Jac(RPP16EE$Com)
\dontrun{
if(require(ade4) && require(adegraphics)){
oldparamadeg <- adegpar()
adegpar("plegend.drawKey" = FALSE)
table.value(RPP16EE$Com)
adegpar(oldparamadeg)
# In this graph, black squares indicate 
# which species (S1-S15) are present in which plot (P1-P9)
}
}
}
\keyword{datasets}
