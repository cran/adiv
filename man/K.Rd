\name{K}
\alias{K}
\alias{Kstar}
\alias{Kw}
\title{
Indices K, Kstar and Kw of Phylogenetic Signal
}
\description{
The function \code{K} calculates Blomberg et al. (2003) statistic \emph{K} and performs a permutation test where species identities are maintained in the phylogeny while the trait values of species are randomly shuffled (permuted) (Pavoine and Ricotta 2013).

The function \code{Kstar} calculates Blomberg et al. (2003) statistic \emph{K*} and performs a permutation test where species identities are maintained in the phylogeny while the trait values of species are randomly shuffled (permuted) (Pavoine and Ricotta 2013).

The function \code{Kw} calculates Pavoine and Ricotta (2013) statistic \eqn{K_w} and performs a permutation test where species identities are maintained in the phylogeny while the trait values of species are randomly shuffled (permuted).
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
a numeric: the number of randomizations
}
  \item{alter}{
a character specifying the alternative hypothesis; must be one of \code{"greater"} (phylogenetic signal; the default), \code{"less"} or \code{"two-sided"}. If several values are given, only the first one is retained.
}
}
\value{
Each function returns an object of class \code{randtest} with the results of the permutation tests. (see function \code{\link{randtest}} in package ade4)
}
\references{
Blomberg, S.P., Garland, T., Ives, A.R. (2003) Testing for phylogenetic signal in comparative data: behavioral traits are more labile. \emph{Evolution}, \bold{57}, 717--745.

Pavoine, S., Ricotta, C. (2013) Testing for phylogenetic signal in biological traits: the ubiquity of cross-product statistics. \emph{Evolution}, \bold{67}, 828--840.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
See Also \code{\link{Kstar}}, \code{\link{Kw}}
}
\examples{
data(rockfish)
phy <- read.tree(text=rockfish$tre)
K(phy, rockfish$traits[phy$tip.label, 1])
Kstar(phy, rockfish$traits[phy$tip.label, 1])
Kw(phy, rockfish$traits[phy$tip.label, 1])
}
\keyword{models}

