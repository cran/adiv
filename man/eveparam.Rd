\name{eveparam}
\alias{eveparam}
\alias{plot.eveparam}
\title{
Parametric Indices of Species Evenness
}
\description{
The function \code{eveparam} calculates parametric evenness indices. The parameter controls the relative importance given to rare versus abundant species in a community.
The function \code{plot.eveparam} plots the results of function \code{eveparam}.
}
\usage{
eveparam(comm, method = c("hill", "tsallis", "renyi"), q = 2, tol = 1e-08)

\method{plot}{eveparam}(x, legend = TRUE, 
legendposi = "topright", axisLABEL = "Evenness", type = "b", 
col = if (is.numeric(x)) NULL 
else sample(colors(distinct = TRUE), nrow(x$eve)), 
lty = if (is.numeric(x)) NULL else rep(1, nrow(x$eve)), 
pch = if (is.numeric(x)) NULL else rep(19, nrow(x$eve)), \dots)
}
\arguments{
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and abundance as entry. 
}
  \item{method}{
a string: either "hill" for the Hill numbers (Hill 1973), "tsallis" for the Tsallis or HCDT entropy (Harvda and Charvat 1967; Daroczy 1970; Tsallis 1988), or "renyi" for Renyi's entropy (Renyi 1960). These indices are divided by the value they would have if all species had even abundances.
}
  \item{q}{
a positive numeric or a vector of positive numerics that give values for the \code{q} parameter.
}
  \item{tol}{
numeric tolerance threshold: values between -tol and tol are considered equal to zero.
}
  \item{x}{
an object of class \code{eveparam} obtained with function \code{eveparam}.
}
  \item{legend}{
a logical. If TRUE a legend is given with the colour, the type of line (etc.) used to define the evenness curve of each community.
}
  \item{legendposi}{
a string or a numeric that gives the position of the legend to be passed to function \code{legend} of the base of R.
}
\item{axisLABEL}{
a string to display on the main axis of the plot to designate what we are measuring. The default is \code{"Evenness"}. 
}
  \item{type}{
a string to be passed to the graphic argument \code{type} of functions \code{plot} and \code{lines} used to draw the evenness curve of each community.
}
  \item{col}{
colours to be passed to the graphic argument \code{col} of functions \code{plot} and \code{lines} to define the colour of the evenness curve of each community.
}
  \item{lty}{
type of line (plain, broken etc.) to be passed to the graphic argument \code{lty} of functions \code{plot} and \code{lines} used to draw the evenness curve of each community.
}
  \item{pch}{
type of point (open circle, close circle, square etc.) to be passed to the graphic argument \code{pch} of functions \code{plot} and \code{lines} used to draw the evenness curve of each community.
}
  \item{\dots}{
other arguments can be added and passed to the functions \code{plot} and \code{lines} used to draw the graphic.
}
}
\value{
Function \code{eveparam}, if only one value of \code{q} is given, returns a vector with the evenness in the communities. If more than one value of \code{q} is given, a list of two objects is returned: 
\item{q}{the vector of values for \code{q};}
\item{div}{a data frame with the evenness in the communities calculated for all values of \code{q}.}

The function \code{plot.eveparam} returns a graphic.
}
\references{
Daroczy, Z. (1970) Generalized information functions. \emph{Information and Control}, \bold{16}, 36--51.

Havrda, M., Charvat, F. (1967) Quantification method of classification processes: concept of structural alpha-entropy. \emph{Kybernatica}, \bold{3}, 30--35.

Hill, M.O. (1973) Diversity and evenness: a unifying notation and its consequences. \emph{Ecology}, \bold{54}, 427--432.

Renyi, A. (1960) On measures of entropy and information. \emph{Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability}, \bold{1}, 547--561.

Tsallis, C. (1988) Possible generalization of Boltzmann-Gibbs statistics. \emph{Journal of Statistical Physics}, \bold{52}, 480--487.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
data(batcomm)
ab <- batcomm$ab
plot(eveparam(ab))
plot(eveparam(ab, q=0:4))
}
\keyword{models}
