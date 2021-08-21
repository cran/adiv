\name{evoeveparam}
\alias{evoeveparam}
\alias{plot.evoeveparam}
\title{
Parametric Indices of Phylogenetic evenness
}
\description{
Function \code{evoeveparam} calculates phylogenetic evenness (evenness in features, which are branch units of a phylogenetic tree) in communities. It uses parametric indices derived from Tsallis and Hill compositional indices, and named qfeveHCDT, qfeveHill, and qfeveRenyi in Pavoine and Ricotta (2019). \code{evoeveparam} can also be applied to functional trees rather than phylogenies, to calculate a functional evenness.
The function \code{plot.evoeveparam} plots the results of function \code{evoeveparam}.
}
\usage{
evoeveparam(phyl, comm, method = c("hill", "tsallis", "renyi"), 
    q = 2, option = 1:3, H = NULL, tol = 1e-8)

\method{plot}{evoeveparam}(x, legend = TRUE, 
    legendposi = "topright", axisLABEL = "Tree-based evenness", 
    type="b", col = if(is.numeric(x)) NULL 
    else sample(colors(distinct = TRUE), nrow(x$eve)), 
    lty = if(is.numeric(x)) NULL else rep(1, nrow(x$eve)), 
    pch = if(is.numeric(x)) NULL else rep(19, nrow(x$eve)), 
    \dots)
}
\arguments{
  \item{phyl}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.}
  \item{comm}{a data frame or a matrix typically with communities as rows, species as columns and abundance as entry. Species should be labeled as in the phylogenetic tree (object \code{phyl}) where they are the tips.}
  \item{method}{a string: either "hill" for qfeveHill using the Hill numbers (Hill 1973), "tsallis" for qfeveHCDT using the Tsallis or HCDT entropy (Harvda and Charvat 1967; Daroczy 1970; Tsallis 1988), or "renyi" for qfeveRenyi using Renyi's entropy (Renyi 1960). If several values are given, only the first one is considered. See details.}
  \item{q}{a vector with nonnegative value(s) for parameter \code{q} of functions qfeveHCDT, qfeveHill, and qfeveRenyi.}
  \item{option}{an integer (either 1, 2 or 3). If 1, the (Hill, Tsallis or Renyi) diversity index is divided by the value it would have if species had same abundance and were independent, if 2, the diversity is divided by the value it would have if species had same abundance, were independent, and at the maximum observed distance from tip to root, if 3, the diversity is divided by the value it would have if species had same abundance, were independent, and at a distance equal to H from the root of the tree. Options 1 and 2 are equivalent in case of an ultrametric tree.}
  \item{H}{a numeric; \code{H} must be higher than the largest observed distance from tip to root.}
  \item{tol}{numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.}
  \item{x}{an object of class \code{evoeveparam} obtained with function \code{evoeveparam}.}
  \item{legend}{a logical. If TRUE a legend is given with the colour, the type of line (etc.) used to define the evenness curve of each community.}
  \item{legendposi}{a string that gives the position of the legend to be passed to function \code{legend} of the base of R.}
\item{axisLABEL}{
a string to display on the main axis of the plot to designate what we are measuring. The default is \code{"Tree-based evenness"}.} 
  \item{type}{a string to be passed to the graphic argument \code{type} of functions \code{plot} and \code{lines} used to draw the evenness curve of each community.}
  \item{col}{vector of colours to be passed to the graphic argument \code{col} of functions \code{plot} and \code{lines} to define the colour of the evenness curve of each community.}
  \item{lty}{vector of type of line (plain, broken etc.) to be passed to the graphic argument \code{lty} of functions \code{plot} and \code{lines} used to draw the evenness curve of each community.}
  \item{pch}{type of point (open circle, close circle, square etc.) to be passed to the graphic argument \code{pch} of functions \code{plot} and \code{lines} used to draw the evenness level of each community.}
  \item{\dots}{other arguments can be added and passed to the functions \code{plot} and \code{lines} used to draw the graphic.}
}
\details{
Function \code{evoeveparam} calculates feature evenness (features = branch units on a phylogenetic [or functional] tree) using parametric indices qfeveHCDT (with \code{method=tsallis}), qfeveHill (with \code{method=hill}), qfeveRenyi (with \code{method=renyi}) developed in Pavoine and Ricotta (2019). Note that Pavoine and Ricotta (2019) recommend the use of index qfeveHill (with \code{method=hill}). 
}
\value{
If only one value of \code{q} is given, the function \code{evoeveparam} returns a vector with the evenness values for the communities.
If more than one value of \code{q} is given, a list of two objects is returned: 
\item{q}{the vector of values for \code{q};}
\item{eve}{a data frame with the phylogenetic evenness in each community calculated for all values of \code{q}.}

The function \code{plot.evoeveparam} returns a graphic.
}
\references{
The methodologies and scripts were developed by

Pavoine, S., Ricotta, C. (2019) A simple translation from indices of species diversity to indices of phylogenetic diversity. \emph{Ecological Indicators}, \bold{101}, 552--561.

using earlier work by:

Chao, A., Chiu, C.-H., Jost, L. (2010) Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society London Series B}, \bold{365}, 3599--3609.

Daroczy, Z. (1970) Generalized information functions. \emph{Information and Control}, \bold{16}, 36--51.

Havrda, M., Charvat F. (1967) Quantification method of classification processes: concept of structural alpha-
entropy. \emph{Kybernetik}, \bold{3}, 30--35.

Hill, M.O. (1973) Diversity and evenness: a unifying notation and its consequences. \emph{Ecology}, \bold{54}, 427--432.

Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 

Renyi, A. (1960) On measures of entropy and information. \emph{Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability}, \bold{1}, 547--561.

Tsallis, C. (1988) Possible generalization of Boltzmann-Gibbs statistics. \emph{Journal of Statistical Physics}, \bold{52}, 480--487.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evodivparam}}, \code{\link{evouniparam}}
}
\examples{
\dontrun{
if(require(ape)){

data(batcomm)
phy <- read.tree(text=batcomm$tre)
ab <- batcomm$ab[, phy$tip.label]
plot(evoeveparam(phy, ab))
plot(evoeveparam(phy, ab, q=seq(0, 10, length=20)))

}
}
}
\keyword{models}
