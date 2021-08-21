\name{evouniparam}
\alias{evouniparam}
\alias{plot.evouniparam}
\title{
Parametric Indices of Phylogenetic Uniqueness
}
\description{
Function \code{evouniparam} calculates phylogenetic uniqueness in communities using parametric indices derived from Tsallis and Hill compositional indices. \code{evouniparam} can also be applied to functional trees rather than phylogenies, to calculate a functional uniqueness.
The function \code{plot.evouniparam} plots the results of function \code{evouniparam}.
}
\usage{
evouniparam(phyl, comm, 
    method = c("hill", "tsallis", "renyi"), 
    q = 2, tol = 1e-08)

\method{plot}{evouniparam}(x, legend = TRUE, 
    legendposi = "topright", 
    axisLABEL = "Tree-based uniqueness",
    type="b", col = if(is.numeric(x)) NULL 
    else sample(colors(distinct = TRUE), nrow(x$uni)), 
    lty = if(is.numeric(x)) NULL else rep(1, nrow(x$uni)), 
    pch = if(is.numeric(x)) NULL else rep(19, nrow(x$uni)), 
    \dots)
}
\arguments{
  \item{phyl}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.}
  \item{comm}{a data frame or a matrix typically with communities as rows, species as columns and abundance as entry. Species should be labeled as in the phylogenetic tree where they are the tips.}
  \item{method}{a string: either "hill" for the Hill numbers (Hill 1973), "tsallis" for the Tsallis or HCDT entropy (Harvda and Charvat 1967; Daroczy 1970; Tsallis 1988), or "renyi" for Renyi's entropy (Renyi 1960). If several values are given, only the first one is considered. See details.}
  \item{q}{a vector with nonnegative value(s) for parameter \code{q}. See details.}
  \item{tol}{numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.}
  \item{x}{an object of class \code{evouniparam} obtained with function \code{evouniparam}.}
  \item{legend}{a logical. If TRUE a legend is given with the colour, the type of line (etc.) used to define the uniqueness curve of each community.}
  \item{legendposi}{a string that gives the position of the legend to be passed to function \code{legend} of the base of R.}
\item{axisLABEL}{
a string to display on the main axis of the plot to designate what we are measuring. The default is \code{"Tree-based uniqueness"}.} 
  \item{type}{a string to be passed to the graphic argument \code{type} of functions \code{plot} and \code{lines} used to draw the uniqueness curve of each community.}
  \item{col}{vector of colours to be passed to the graphic argument \code{col} of functions \code{plot} and \code{lines} to define the colour of the uniqueness curve of each community.}
  \item{lty}{vector of type of line (plain, broken etc.) to be passed to the graphic argument \code{lty} of functions \code{plot} and \code{lines} used to draw the uniqueness curve of each community.}
  \item{pch}{type of point (open circle, close circle, square etc.) to be passed to the graphic argument \code{pch} of functions \code{plot} and \code{lines} used to draw the uniqueness level of each community.}
  \item{\dots}{other arguments can be added and passed to the functions \code{plot} and \code{lines} used to draw the graphic.}
}
\details{
Function \code{evouniparam} calculates feature uniqueness (features = branch units on a phylogenetic [or functional] tree) using parametric indices qfuniHCDT (with \code{method=tsallis}), qfuniHill (with \code{method=hill}), qfuniRenyi (with \code{method=renyi}) developed in Pavoine and Ricotta (2019). Note that Pavoine and Ricotta (2019) recommend the use of index qfuniHill (with \code{method=hill}). 
}
\value{
If only one value of \code{q} is given, the function \code{evouniparam}  returns a vector with the phylogenetic uniqueness of each community.
If more than one value of \code{q} is given, a list of two objects is returned: 
\item{q}{the vector of values for \code{q};}
\item{uni}{a data frame with the phylogenetic uniqueness in each community calculated for all values of \code{q}.}

The function \code{plot.evouniparam} returns a graphic.
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
\code{\link{evodivparam}}, \code{\link{evoeveparam}}
}
\examples{
\dontrun{
if(require(ape)){

data(batcomm)
phy <- read.tree(text=batcomm$tre)
ab <- batcomm$ab[, phy$tip.label]
plot(evouniparam(phy, ab))
plot(evouniparam(phy, ab, q=seq(0, 10, length=20)))

}
}
}
\keyword{models}
