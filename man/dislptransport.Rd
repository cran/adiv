\name{dislptransport}
\alias{dislptransport}
\title{Algorithmic index of plot-to-plot functional (or phylogenetic) dissimilarity}
\description{
The function \code{dislptransport} calculates the DKG measure of dissimilarity (Kosman 1996; Gregorius et al. 2003) applied here to functional differences between plots (eq. 3 in Ricotta et al. 2021).
}
\usage{
dislptransport(comm, dis, diag = FALSE, upper = FALSE, Nind = 10000)
}
\arguments{
\item{comm}{a matrix containing the relative or absolute abundance of all species in plots. Columns are species and plots are rows. Column labels (species names) should be assigned as in \code{dis}.}
\item{dis}{a matrix or an object of class dist providing the functional dissimilarities between species (dissimilarities are nonnegative, symmetric, and the dissimilarity between a species and itself is zero). \bold{Species here must be in the same order as in the columns of comm}.}
\item{diag}{logical value indicating whether the diagonal of the dissimilarity matrix should be printed by \code{print.dist}.}
\item{upper}{logical value indicating whether the upper triangle of the distance matrix should be printed by \code{print.dist}.}
\item{Nind}{an integer. The algorithmic index will be applied by assuming that each plot contains \code{Nind} individuals. The highest \code{Nind}, the most precise the index value will be (see Gregorius et al. 2003, for more details).}
}
\value{
The function \code{dislptransport} returns an object of class \code{dist} matrix with the values of the dissimilarity index DKG (Ricotta et al. 2021) for each pair of plots. 
}
\references{
Ricotta, C., Kosman, E., Laroche, F., Pavoine, S. (2021) Beta redundancy for functional ecology. \emph{Methods in Ecology and Evolution}, \bold{12}, 1062--1069. \doi{10.1111/2041-210X.13587}

Gregorius, H.-R., Gillet, E.M., Ziehe, M. (2003) Measuring differences of trait distributions between populations. \emph{Biometrical Journal}, \bold{8}, 959--973. \doi{10.1002/bimj.200390063}

Kosman, E. (1996) Difference and diversity of plant pathogen populations: a new approach for measuring. \emph{Phytopathology}, \bold{86}, 1152--1155.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{betaUniqueness}} for plot-to-plot functional (or phylogenetic) uniqueness indices.
}
\examples{
\dontrun{
data(RutorGlacier)
fundis <- dist(scale(RutorGlacier$Traits2[1:6]))
fundis <- fundis/max(fundis)
DKG <- dislptransport(RutorGlacier$Abund, fundis)
}
}
\keyword{models}

