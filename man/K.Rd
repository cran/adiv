\name{K}
\alias{K}
\alias{Kstar}
\alias{Kw}
\title{
Indices K, Kstar and Kw of Phylogenetic Signal
}
\description{
The functions \code{K}, \code{Kstar} and \code{Kw} calculate Blomberg et al. (2003) statistics \emph{K}, and \emph{K*} and Pavoine and Ricotta (2013) statistic \eqn{K_w}{Kw}, respectively. Then they perform a permutation test where species identities are maintained in the phylogeny while the trait values of species are randomly shuffled (permuted) (Pavoine and Ricotta 2013).
}
\usage{
K(phyl, trait, nrep = 999, alter = c("greater", "less", "two-sided"))

Kstar(phyl, trait, nrep = 999, alter = c("greater", "less", "two-sided"))

Kw(phyl, trait, nrep = 999, alter = c("greater", "less", "two-sided"))
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}.
}
  \item{trait}{
a vector with the trait value for each species (tip) in the phylogenetic tree. Trait values for species must be in the same order as species in the phylogenetic tree.
}
  \item{nrep}{
a numeric: the number of randomizations.
}
  \item{alter}{
a string specifying the alternative hypothesis; it must be one of \code{"greater"} (phylogenetic signal; the default), \code{"less"} or \code{"two-sided"}. If several values are given, only the first one is retained.
}
}
\value{
Each function returns an object of class \code{randtest} with the results of the permutation tests. (see function \code{randtest} in package ade4)
}
\details{
Blomberg et al. (2003) introduced two statistics of phylogenetic signal:
\deqn{K = MSE0 / MSE}
\deqn{K^* = MSE^* / MSE}{K* = MSE* / MSE}
where \emph{MSE} is the mean squared error of the trait values calculated using the variance-covariance matrix derived from the phylogenetic tree, \emph{MSE0} is the mean squared error of the tip trait values, measured from a phylogenetically correct mean of tip trait values and \emph{MSE*} is the mean squared error of the tip trait values, measured from the estimate of the mean of the raw tip trait values. In both statistics \emph{K} and \emph{K*}, the value of \emph{MSE} will be relatively small if the phylogenetic tree accurately describes the variance-covariance pattern observed in the data, leading to high values for \emph{K} and \emph{K*} (meaning high phylogenetic signal). In functions \code{K} and \code{Kstar}, \emph{K} and \emph{K*}, respectively, are divided (normalized) by their expected value if the trait evolved under a Brownian motion along the branches of the phylogenetic tree (this expected value is invariant under permutation of trait values among the tips of the phylogeny).

To test for phylogenetic signal, Blomberg et al. (2003) actually considered neither \emph{K} nor \emph{K*} but \emph{MSE} as the core statistic associated with random permutations of trait values among tips of the phylogenetic tree. Although the literature on phylogenetic signal has currently mostly ignored \emph{K*} focusing on statistic \emph{K}, \emph{K*} could thus actually have been considered as the core statistic of Blomberg et al. (2003) test for phylogenetic signal. Indeed, as \emph{MSE*} is independent of permutations of trait values among the tips of the phylogeny while \emph{MSE0} is, Blomberg et al. (2003) approach corresponds to considering \emph{K*} and not \emph{K} as the statistic of the test of phylogenetic signal in traits. This test is also equivalent to an alternative implemented via phylogenetically independent contrasts also proposed by Blomberg et al. (2003).

Function \code{KW} implements index \eqn{K_w}{Kw}, a modified version of \emph{K*} that grants a higher importance in the calculation of phylogenetic signal to the tips that have many closely related tips (Pavoine and Ricotta 2013).

In functions, \code{K}, \code{Kstar} and \code{Kw}, I considered the same permutation scheme as in Blomberg et al. (2003) but used \emph{K}, \emph{K*} and \eqn{K_w}{Kw}, as the core statistic, respectively. The test developed by Blomberg et al. (2003) thus corresponds to function \code{Kstar}.
}
\references{
Blomberg, S.P., Garland, T., Ives, A.R. (2003) Testing for phylogenetic signal in comparative data: behavioral traits are more labile. \emph{Evolution}, \bold{57}, 717--745.

Pavoine, S., Ricotta, C. (2013) Testing for phylogenetic signal in biological traits: the ubiquity of cross-product statistics. \emph{Evolution}, \bold{67}, 828--840.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{Kstar}}, \code{\link{Kw}}
}
\examples{
\dontrun{
if(require(ape) && require(ade4)){

data(rockfish)
phy <- read.tree(text=rockfish$tre)

theK <- K(phy, rockfish$traits[phy$tip.label, 1])
theK
plot(theK)

theKstar <- Kstar(phy, rockfish$traits[phy$tip.label, 1])
theKstar
plot(theKstar)

theKw <- Kw(phy, rockfish$traits[phy$tip.label, 1])
theKw
plot(theKw)

}
}
}
\keyword{models}

