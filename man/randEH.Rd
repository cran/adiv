\name{randEH}
\alias{randEH}
\title{
Amount of Evolutionary History Preserved when Random Species are Saved
}
\description{
When branch lengths in an ultrametric phylogenetic tree are expressed as divergence times, the total sum of branch lengths in that tree expresses the amount of evolutionary history. The function \code{randEH} calculates the amount of evolutionary history preserved when \emph{k} random species out of \emph{n} original species are saved. 
}
\usage{
randEH(phyl, nbofsp, nrep = 10)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}.
}
  \item{nbofsp}{
an integer indicating the number of species saved (\emph{k}).}
  \item{nrep}{
an integer indicating the number of random sampling.}
}
\value{
Function \code{randEH} returns a numeric vector with the amount of evolutionary history preserved by each random drawing of the \emph{k} species to be saved.
}
\references{
Nee, S. and May, R.M. (1997) Extinction and the loss of evolutionary history. \emph{Science}, \bold{278}, 692--694.

Pavoine, S., Ollier, S. and Dufour, A.-B. (2005) Is the originality of a species measurable? \emph{Ecology Letters}, \bold{8}, 579--586.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
with contributions of Stephane Dray
}
\seealso{
\code{\link{optimEH}}
}
\examples{
\dontrun{
if(require(ape) && require(adephylo)){

data(carni70, package = "adephylo")
tre <- read.tree(text=carni70$tre)
adiv:::EH(tre)
R <- adiv:::randEH(tre, 10, nrep=1000)
hist(R)

}
}
}
\keyword{models}
