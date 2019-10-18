\name{optimEH}
\alias{optimEH}
\title{
Nee and May's Optimizing Process 
}
\description{
The function \code{optimEH} performs Nee and May's optimizing scheme. When branch lengths in an ultrametric phylogenetic tree are expressed as divergence times, the total sum of branch lengths in that tree expresses the amount of evolutionary history. Nee and May's algorithm optimizes the amount of evolutionary history preserved if only \emph{k} species out of \emph{n} were to be saved. The \emph{k}-1 closest-to-root nodes are selected, which defines \emph{k} clades; one species from each clade is picked. At this last step, we decide to select the most original species of each of the \emph{k} clades. }
\usage{
optimEH(phyl, nbofsp, tol = 1e-08, give.list = TRUE)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}.
}
  \item{nbofsp}{
an integer indicating the number of species saved (\emph{k}).
}
  \item{tol}{
a tolerance threshold for null values (a value less than \code{tol} in absolute terms is considered as zero). }
  \item{give.list}{
a logical indicating whether a list of optimizing species should be provided. If \code{give.list = TRUE}, \code{optimEH} provides the list of the \emph{k} species that optimize the amount of evolutionary history preserved and are the most original species in their clades. If \code{give.list = FALSE}, \code{optimEH} returns directly the real value giving the amount of evolutionary history preserved.}
}
\value{
If \code{give.list=FALSE}, function \code{optimEH} returns a numeric.

If \code{give.list=TRUE}, function \code{optimEH} returns a list containing:
\item{value}{a real value providing the amount of evolutionary history preserved;}
\item{selected.sp}{a data frame containing the list of the k species which optimize the amount of evolutionary history preserved and are the most original species in their clades.}
}
\references{
Pavoine, S., Ollier, S. and Dufour, A.-B. (2005) Is the originality of a species measurable? \emph{Ecology Letters}, \bold{8}, 579--586. 
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
with contributions of Stephane Dray
}
\seealso{
\code{\link{randEH}}
}
\examples{
\dontrun{
if(require(ape) && require(adephylo)){

data(carni70, package="adephylo")
tre <- read.tree(text=carni70$tre)
adiv:::EH(tre)
adiv:::optimEH(tre, 10)

}
}
}
\keyword{models}
