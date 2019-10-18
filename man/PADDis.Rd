\name{PADDis}
\alias{PADDis}
\alias{DJac}
\alias{Jac}
\title{
Functional Dissimilarity Measures for Presence-Absence Data}
\description{
Functions \code{PADDis}, \code{DJac} and \code{Jac} calculate the dissimilarity coefficients introduced in Ricotta et al. (2016). These dissimilarity coefficients use traditional mismatching components a, b and c of the 2 x 2 contingency table expressed as to include functional or phylogenetic differences among species and noted A, B, C. Components B and C represent the functional or phylogenetic uniqueness of community X compared with community Y and vice versa. Component A represents the functional or phylogenetic similarities between communities X and Y.}
\usage{
PADDis(comm, dis, method = NULL, diag = FALSE, upper = FALSE)

DJac(comm, dis, diag = FALSE, upper = FALSE)

Jac(comm, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{comm}{
a matrix or a data frame with communities (or plots, assemblages, etc.) as rows and species as columns containing the incidence (0/1) of all species in the communities.}
  \item{dis}{
an object of class \code{dist} containing the (functional) dissimilarities among species.}
  \item{method}{
an integer between 0 and 6. If \code{NULL} the choice is made with a console message. See details.}
  \item{diag}{
a logical value indicating whether the diagonal of the distance matrix should be printed by function \code{print.dist}.}
  \item{upper}{
a logical value indicating whether the upper triangle of the distance matrix should be printed by function \code{print.dist}.}
}
\details{
In \code{PADDIS}, dissimilarities among communities are calculated with the following formulas:

Generalized Jaccard dissimilarity, with method = 1
\deqn{\frac{B+C}{a+b+c}}{(B+C)/(a+b+c)}

Generalized Sorensen dissimilarity, with method = 2
\deqn{\frac{B+C}{2a+b+c}}{(B+C)/(2a+b+c)}

Generalized Sokal and Sneath dissimilarity, with method = 3
\deqn{\frac{2(B+C)}{a+2(b+c)}}{[2(B+C)]/[a+2(b+c)]}

Generalized Ochiai dissimilarity, with method = 4
\deqn{\frac{\sqrt{A+B}\sqrt{A+C}-A}{\sqrt{a+b}\sqrt{a+c}}}{[sqrt(A+B)sqrt(A+C)-A]/[sqrt(a+b)sqrt(a+c)]}

Generalized Simpson dissimilarity, with method = 5
\deqn{\frac{min\{B+C\}}{a+min\{b+c\}}}{[min(B+C)][a+min(b+c)]}

Generalized Kulczynski dissimilarity, with method = 6
\deqn{0.5*(\frac{B}{a+b}+\frac{C}{a+c})}{0.5*(B/(a+b)+C/(a+c))}

\code{DJac} and \code{Jac} use the additive decomposition of the Jaccard index into turnover and richness difference. \code{DJac} takes into account the (functional or phylogenetic) dissimilarities among species while \code{Jac} does not.
}
\value{
In function \code{PADDis}, if \code{method=0}, then the function \code{PADDis} returns 6 matrices corresponding to the a, b, c, A, B, and C values per pair of communities. Otherwise, it returns an object of class \code{dist} corresponding to the dissimilarities among communities.

Functions \code{DJac} and \code{Jac} return a list of three objects of class \code{dist}:
\item{J}{for the full dissimilarities between communities;}
\item{JRepl}{for the turnover component of the dissimilarities;}
\item{JRich}{for the component of difference in richness.}
}
\references{
Ricotta, C., Podani, J., Pavoine, S. (2016) A family of functional dissimilarity measures for presence and absence data. \emph{Ecology and Evolution}, \bold{6}, 5383--5389}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
data(RPP16EE)
Com <- RPP16EE$Com
Dis <- as.dist(RPP16EE$Dis)
J <- Jac(Com)
DJ <- DJac(Com, Dis)

plot(c(as.matrix(DJ$J)[1,]), ylab="Dissimilarity", 
xlab="Plot-to-plot comparison", pch=15, type="b", 
ylim=c(0,1), main="Jaccard")

lines(c(as.matrix(J$J)[1,]), type="b", pch=18)

legend("bottomright", legend=c("P/A scores", "functional data"), 
pch=c(15,18), lty=1)

plot(c(as.matrix(DJ$JRepl)[1,]), ylab="Dissimilarity",
xlab="Plot-to-plot comparison", pch=15, type="b", 
ylim=c(0,1), main="Species replacement")

lines(c(as.matrix(J$JRepl)[1,]), type="b", pch=18)

legend("bottomright", legend=c("P/A scores", "functional data"), 
pch=c(15,18), lty=1)


#Use the following instruction to obtain all components:

PADDis(Com, Dis, method=0)
}
\keyword{models}
