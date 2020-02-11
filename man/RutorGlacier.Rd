\name{RutorGlacier}
\alias{RutorGlacier}
\docType{data}
\title{
Functional and phylogenetic composition of plant communities along a primary succession on glacial deposits. 
}
\description{
59 plots, each of about 25 m2 in size, were sampled above the tree line at the foreland of the Rutor Glacier (Italy). The abundance of plant species were collected by Caccianiga et al. (2006); Ricotta et al. (2016) established trait data using Grime's (1974) plant strategy theory and Ricotta et al. (2015) the phylogenetic tree using the dated Daphne phylogeny (Durka and Michalski 2012) for European flora.
}
\usage{data("RutorGlacier")}
\format{
\code{RutorGlacier} is a list of four components:

\code{Abund}, a data frame with plots as rows, species as columns and abundance of species in plots as entries.

\code{Traits}, a data frame with species as rows and three traits (named C, S, and R) as columns. Species were classified in terms of Grime's (1974) plant strategy theory, as competitors (C), stress tolerators (S) and ruderals (R) with fuzzy-coded values in the range 0-100, such that the sum of C+S+R was equal to 100.

\code{TreeNW}, a phylogenetic tree in newick format for all plant species in \code{Abund} and \code{traits}.

\code{Fac}, a vector of strings defining in which successional stage each plot belongs to: "early" = early-successional stage, "mid" = mid-successional stage and "late" = late-successional stage.
}
\references{
Caccianiga, M., Luzzaro, A., Pierce, S., Ceriani, R.M., Cerabolini, B.E.L. (2006) The functional basis of a primary succession resolved by CSR classification. \emph{Oikos}, \bold{112}, 10--20.

Durka, W., Michalski, S.G. (2012) Daphne: a dated phylogeny of a large European flora for phylogenetically informed ecological analyses. \emph{Ecology}, \bold{93}, 2297.

Grime, J.P. (1974) Vegetation classification by reference to strategies. \emph{Nature}, \bold{250}, 26--31.

Ricotta, C., Bacaro, G., Caccianiga, M., Cerabolini, B.E.L., Moretti, M. (2015) A classical measure of phylogenetic dissimilarity and its relationship with beta diversity. \emph{Basic and Applied Ecology}, \bold{16}, 10--18.

Ricotta, C., de Bello, F., Moretti, M., Caccianiga, M., Cerabolini, B.E.L., Pavoine, S. (2016) Measuring the functional redundancy of biological communities: a quantitative guide. \emph{Methods in Ecology and Evolution}, \bold{7}, 1386--1395.

Ricotta, C., Laroche, F., Szeidl, L., Pavoine, S. (2020) From alpha to beta functional and phylogenetic redundancy. \emph{Methods in Ecology and Evolution}. In press.
}
\examples{
\dontrun{
if(require(ape) && require(ade4)){
data(RutorGlacier)
phy <- read.tree(text=RutorGlacier$TreeNW)
plot(phy)
ab <- RutorGlacier$Abund[, phy$tip.label]
plot(abgevodivparam(phy, ab, q=0:4))

# Phylogenetic dissimilarities between plots
# (Ricotta et al. 2020)
Dp <- DP(phy, ab, tol=0.00001)
pcoDp <- dudi.pco(sqrt(Dp), full=TRUE)
s.class(pcoDp$li, as.factor(RutorGlacier$Fac))

# Phylogenetic beta uniqueness (Ricotta et al. 2020)
Up <- betaTreeUniqueness(phy, ab, tol=0.00001)
# Average uniqueness between two plots at each successional stage
fac <- as.factor(RutorGlacier$Fac)
mean(Up[fac == "early", fac == "early"])
mean(Up[fac == "mid", fac == "mid"])
mean(Up[fac == "late", fac == "late"])
}
}
}
\keyword{datasets}
