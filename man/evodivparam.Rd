\name{evodivparam}
\alias{evodivparam}
\alias{plot.evodivparam}
\title{
Parametric Indices of Phylogenetic Diversity
}
\description{
Function \code{evodivparam} calculates phylogenetic diversity in communities using parametric indices derived from Tsallis and Hill compositional indices. It can also be applied to functional trees rather than phylogenies, to calculate a functional diversity.
The function \code{plot.evodivparam} plots the results of function \code{evodivparam}.
}
\usage{
evodivparam(phyl, comm, 
    method = c("hill", "tsallis", "renyi"), 
    q = 2, tol = 1e-08)

\method{plot}{evodivparam}(x, legend = TRUE, 
    legendposi = "topright", axisLABEL = "Tree-based diversity", 
    type="b", col = if(is.numeric(x)) NULL 
    else sample(colors(distinct = TRUE), nrow(x$div)), 
    lty = if(is.numeric(x)) NULL else rep(1, nrow(x$div)), 
    pch = if(is.numeric(x)) NULL else rep(19, nrow(x$div)), 
    \dots)
}
\arguments{
  \item{phyl}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.}
  \item{comm}{a data frame or a matrix typically with communities as rows, species as columns and abundance as entry. Species should be labeled as in the phylogenetic tree where they are the tips.}
  \item{method}{a string: either "hill" for the Hill numbers (Hill 1973), "tsallis" for the Tsallis or HCDT entropy (Harvda and Charvat 1967; Daroczy 1970; Tsallis 1988), or "renyi" for Renyi's entropy (Renyi 1960). If several values are given, only the first one is considered. See details.}
  \item{q}{a vector with nonnegative value(s) for parameter \code{q}. See details.}
  \item{tol}{numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.}
  \item{x}{an object of class \code{evodivparam} obtained with function \code{evodivparam}.}
  \item{legend}{a logical. If TRUE a legend is given with the colour, the type of line (etc.) used to define the diversity curve of each community.}
  \item{legendposi}{a string that gives the position of the legend to be passed to function \code{legend} of the base of R.}
\item{axisLABEL}{
a string to display on the main axis of the plot to designate what we are measuring. The default is \code{"Tree-based diversity"}.} 
  \item{type}{a string to be passed to the graphic argument \code{type} of functions \code{plot} and \code{lines} used to draw the diversity curve of each community.}
  \item{col}{vector of colours to be passed to the graphic argument \code{col} of functions \code{plot} and \code{lines} to define the colour of the diversity curve of each community.}
  \item{lty}{vector of type of line (plain, broken etc.) to be passed to the graphic argument \code{lty} of functions \code{plot} and \code{lines} used to draw the diversity curve of each community.}
  \item{pch}{type of point (open circle, close circle, square etc.) to be passed to the graphic argument \code{pch} of functions \code{plot} and \code{lines} used to draw the diversity level of each community.}
  \item{\dots}{other arguments can be added and passed to the functions \code{plot} and \code{lines} used to draw the graphic.}
}
\details{
Consider a phylogenetic tree \emph{T}, \eqn{b_T}{b_T} the set of branches in \emph{T}, \emph{k} a branch, \eqn{L_k}{L_k} the length of branch \emph{k}, \emph{j} a community (\emph{j}=1,...,\emph{m}), \eqn{a_{jk}}{a_jk} the abundance associated with branch \emph{k} in community \emph{j} (sum of abundance of all species descending from the branch). \emph{q} is the parameter that increases with the importance given to abundant species compared to rare species in diversity.

The methods available to calculate the phylogenetic diversity in community \emph{j} are:
\code{tsallis}:

\deqn{^q{evoTsallis}_j=
\left[1-\sum_{k \in b_T} L_k \left(\frac{a_{jk}}{\sum_{k \in b_T} L_k a_{jk}}\right)^q \right]/(q-1)}{qevoTsallis=[1-sum_k L_k (a_jk/(sum_k L_k a_jk))^q]/(q-1)}

\code{hill}:

\deqn{^q{evoHill}_j=\left[\sum_{k \in b_T} L_k \left(\frac{a_{jk}}{\sum_{k \in b_T} L_k a_{jk}}\right)^q \right]^{1/(1-q)}}{qevoHill=[sum_k L_k (a_jk/(sum_k L_k a_jk))^q]^(1/(1-q))}

\code{renyi}:
\deqn{^q{evoRenyi}_j=log(^q{evoHill}_j)}{qevoRenyi=log(qevoHill)}

}
\value{
If only one value of \code{q} is given, a vector with the phylogenetic diversity of each community is returned.
If more than one value of \code{q} is given, a list of two objects is returned: 
\item{q}{the vector of values for \code{q};}
\item{div}{a data frame with the phylogenetic diversity of each community calculated for all values of \code{q}.}

The function \code{plot.evodivparam} returns a graphic.
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
\code{\link{abgevodivparam}}, \code{\link{divparam}}
}
\examples{
\dontrun{
if(require(ape)){

data(batcomm)
phy <- read.tree(text=batcomm$tre)
ab <- batcomm$ab[, phy$tip.label]
plot(evodivparam(phy, ab))
plot(evodivparam(phy, ab, q=seq(0, 10, length=20)))

}
}
}
\keyword{models}
