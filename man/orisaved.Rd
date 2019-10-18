\name{orisaved}
\alias{orisaved}
\title{
Maximal or Minimal Amount of Originality Saved under Optimal Conditions
}
\description{
Function \code{orisaved} computes the maximal or minimal amount of originality saved over all combinations of species optimizing the amount of evolutionary history preserved. The originality of a species is measured with the QE-based index (Pavoine et al. 2005).}
\usage{
orisaved(phyl, rate = 0.1, method = 1:2)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}.
}
  \item{rate}{
a real value (between 0 and 1) indicating how many species will be saved for each calculation. For example, if the total number of species is 70 and 'rate = 0.1' then the calculations will be done at a rate of 10\% i.e. for 0 (= 0 \%), 7 (= 10 \%), 14 (= 20 \%), 21 (= 30 \%), ..., 63 (= 90 \%) and 70(= 100 \%) species saved. If 'rate = 0.5' then the calculations will be done for only 0 (= 0 \%), 35 (= 50 \%) and 70(= 100 \%) species saved.
}
  \item{method}{
an integer either 1 or 2 (see details).
}
}
\details{
1 = maximum amount of originality saved 2 = minimum amount of originality saved
}
\value{
Returns a numeric vector. 
}
\references{
Pavoine, S., Ollier, S. and Dufour, A.-B. (2005) Is the originality of a species measurable? \emph{Ecology Letters}, \bold{8}, 579--586. 
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
with contributions of Stephane Dray
}
\seealso{
\code{\link{optimEH}}, \code{\link{randEH}}
}
\examples{
\dontrun{
if(require(ape) && require(adephylo)){

data(carni70, package="adephylo")
tre <- read.tree(text=carni70$tre)
tmax <- adiv:::orisaved(tre, rate = 1 / 70, method = 1) 
tmin <- adiv:::orisaved(tre, rate = 1 / 70, method = 2) 
plot(c(0, 1:70), tmax, xlab = "nb of species saved", ylab = "Originality saved", type = "l") 
lines(c(0, 1:70), tmin, lty = 2) 

}
}
}
\keyword{models}
