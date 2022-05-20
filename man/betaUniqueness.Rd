\name{betaUniqueness}
\alias{betaUniqueness}
\title{Algorithmic index of plot-to-plot functional (or phylogenetic) dissimilarity and uniqueness}
\description{
The function \code{betaUniqueness} calculates uniqueness and redundancy taking account of functional dissimilarities between species using equation 5 and 6 in Ricotta et al. (2021). Note that functional dissimilarities could be replaced by any other type of dissimilarities between species, including phylogenetic dissimilarities.
}
\usage{
betaUniqueness(comm, dis, Nind = 10000)
}
\arguments{
\item{comm}{a matrix containing the relative or absolute abundance of all species in plots. Columns are species and plots are rows. Column labels (species names) should be assigned as in \code{dis}.}
\item{dis}{a matrix or an object of class dist providing the functional dissimilarities between species (dissimilarities are nonnegative, symmetric, and the dissimilarity between a species and itself is zero). \bold{Species here must be in the same order as in the columns of comm}.}
\item{Nind}{an integer. The algorithmic index will be applied by assuming that each plot contains \code{Nind} individuals. The highest \code{Nind}, the most precise the index value will be (see Gregorius et al. 2003, for more details).}
}
\value{
The function \code{betaUniqueness} returns a list with the following objects:

- betaUniqueness: a matrix with the values of the proposed beta uniqueness (\eqn{U\beta}{Ubeta}=DKG/DR) for each pair of plots (Ricotta et al. (2021), eq. 6);

- betaRedundancy: a matrix with the values of the proposed beta redundancy (\eqn{R\beta}{Rbeta}=1-DKG/DR) for each pair of plots (Ricotta et al. (2021), eq. 5);

- dissimilarityGap: a matrix with the values of the dissimilarity gap index (DR-DKG) for each pair of plots;

- DR: a matrix with the values of the species-based (Rogers) dissimilarity index (DR) for each pair of plots (Ricotta et al. (2021), eq. 4);

- DKG: a matrix with the values of the algorithmic functional dissimilarity index (DKG) for each pair of plots (Ricotta et al. (2021), eq. 3).
}
\references{
Ricotta, C., Kosman, E., Laroche, F., Pavoine, S. (2021) Beta redundancy for functional ecology. \emph{Methods in Ecology and Evolution}, \bold{12}, 1062--1069. \doi{10.1111/2041-210X.13587}

Gregorius, H.-R., Gillet, E.M., Ziehe, M. (2003) Measuring differences of trait distributions between populations. \emph{Biometrical Journal}, \bold{8}, 959--973. \doi{10.1002/bimj.200390063}
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{betaTreeUniqueness}} adapted to the use of phylogenetic trees with species as tips, \code{\link{dislptransport}} for the algorithmic functional dissimilarity index (DKG in Ricotta et al. 2021), and \code{\link{uniqueness}} for alpha uniqueness
}
\examples{
\dontrun{
data(RutorGlacier)
fundis <- dist(scale(RutorGlacier$Traits2[1:6]))
fundis <- fundis/max(fundis)
frameDKG <- betaUniqueness(RutorGlacier$Abund, fundis)

f1 <- unlist(sapply(1:58, function(i) rep(RutorGlacier$Fac[i], 59-i)))
f2 <- unlist(sapply(1:58, function(i) RutorGlacier$Fac[-(1:i)]))
f <- paste(f1, f2, sep="-")
F <- factor(f, levels=c("early-early", "mid-mid", "late-late", "early-mid", 
    "mid-late", "early-late"))

vbetaU_A <- as.vector(as.dist(frameDKG$betaUniqueness))

boxplot(vbetaU_A~F, ylab="Beta uniqueness", xlab="Compared successional stages")
# See Ricotta et al. 2021 Electronic Appendix 3 for for details
}
}
\keyword{models}

