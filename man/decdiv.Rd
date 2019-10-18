\name{decdiv}
\alias{decdiv}
\alias{plot.decdiv}
\alias{rtestdecdiv}
\title{
Decomposition of trait-based diversity along the nodes of a phylogenetic tree
}
\description{
The function \code{decdiv} calculates trait-based differences between the lineages that descend from a node of a phylogenetic tree in one or several communities (using presence/absence or abundance data).  

The function \code{plot.decdiv} plots the result of function decdiv for one of the communities. 

The function \code{rtestdecdiv} tests, for one community (with presence/absence or abundance data), if the representation of trait diversity on the phylogenetic tree highlights a nonrandom pattern. 
}
\usage{
decdiv(phyl, comm, dis = NULL, tol = 1e-08, option = 1:5, 
     formula = c("QE", "EDI"))

\method{plot}{decdiv}(x, ncom = 1, col = "black", 
    csize = 1, legend = TRUE, \dots)

rtestdecdiv(phyl, vecab, dis = NULL, tol = 1e-08, 
    option = 1:5, formula = c("QE", "EDI"), 
    vranking = c("complexity", "droot"), 
    ties.method = c("average", "first", "last", "random", 
    "max", "min"), statistic = 1:3, optiontest = NULL, nrep = 99)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust} with species as tips. To ease the interpretation of the results, I advise you to add labels to the nodes of the phylogeny. For example, If your phylogeny is named 'tree' and is of class \code{phylo}, you can name the nodes with the following command:
\code{tree$node.label <- paste("n", 1:tree$Nnode, sep="")}. Use \code{plot(tree, show.node.label=TRUE)} to see the result. If it is of class \code{phylo4}, you can use: \code{nodeLabels(tree) <- paste("n", 1:nNodes(tree), sep="")}. 
}
  \item{comm}{
a vector with species presence/absence or an index of abundance as entries, or a data frame or a matrix typically with communities as rows, species as columns and presence/absence or an index of abundance as entries. Species names in \code{comm} must be similar as in \code{phyl}. Even if \code{comm} is a vector, each entry needs to have a species name (use function \code{name} to attribute names to a vector).
}
  \item{dis}{
either \code{NULL} or an object of class \code{dist} that contains the trait-based distances among species. If \code{NULL}, the Gini-Simpson index is used.
}
  \item{tol}{
a tolerance threshold (a value between -\code{tol} and \code{tol} is considered equal to zero)}
  \item{option}{a numeric (either 1, 2, 3, 4 or 5) indicating the option to use to calculate the contribution of each node of the phylogenetic tree to trait-based diversity. See details. If several values are given, the function retains only the first one.}
  \item{formula}{
either \code{"QE"} (default) or \code{"EDI"}. See details. If several values are given, the function retains only the first one.
}
  \item{x}{
an object of class \code{decdiv}.
}
  \item{ncom}{
if \code{comm} is a matrix, the number of its row corresponding to the focal community for which the results of decdiv will be plotted. Igored if \code{comm} is a vector. 
}
   \item{col}{
the color of circles displayed at each node.
}
   \item{csize}{
a positive numeric giving the scale for plotting the circle at each node. 1 is the default size; if zero, no circle is drawn.
}
   \item{legend}{
a logical indicating whether the legend for the circle size needs to be displayed.
}
  \item{\dots}{
further arguments that can be specified to the internal use of function \code{plot.phylo} (argument \code{y.lim} cannot be modified).
}
   \item{vecab}{a numeric vector giving the presence/absence(1/0) or abundance(non-negative value) of species in a community.}
   \item{vranking}{a string with 2 possible values: either "complexity" for a ranking according to the complexity of the subtree rooted on each interior node (see Pavoine et al. 2010 for an explanation), or "droot" for ordering interior nodes by the distance between them and the root node of the tree.}
   \item{ties.method}{a string to be passed to function \code{rank} of the base of R. It can be one of "average", "first", "last", "random", "max", "min".}
   \item{statistic}{a numeric value or a vector of numeric values. Possible values are 1, 2, or 3. They correspond to the three statistics S1, S2 and S3, respectively, developed by Pavoine et al. (eqs. 5 to 7 in Pavoine et al. 2010).}
   \item{optiontest}{a vector of strings specifying the alternative hypothesis of each test, which must be one of "greater", "less" or "two-sided". If null, then statistic=1 is associated with "greater" and statistic=2 and =3 with "two-sided". See function \code{as.randtest} of package ade4 for details on these alternatives. The length of the vector given to \code{optiontest} must be equal to that given to  \code{statistic}.}
   \item{nrep}{numeric; the number of permutations to be done in each permutation test.}
}
\value{
Function \code{decdiv} returns a matrix with nodes of the phylogenetic tree as rows and the decomposition of trait-based diversity in communities as columns. 
An attribute of this matrix is the phylogenetic tree (of class \code{phylo} with specified names for nodes). If the nodes of \code{phyl} had no names, the function \code{decdiv} automatically attributed names to them. 
}
\details{
The function \code{decdiv} relies on Rao's (1982) quadratic entropy (QE) to measure the trait-based diversity of a set of species. Two formulas for QE have been introduced in the literature one is the original formula by Rao (1982) (which corresponds to \code{formula = "QE"}) and the other one introduced by Champely and Chessel (2002), named Euclidean Diversity Index (which corresponds to \code{formula = "EDI"}). See function \code{\link{QE}} for more details.

In function \code{decdiv}, each node has a weight proportional to the summed relative abundance of its descending species (or to the relative number of descending species if presence/absence data are used). 

With \code{option = 1}, the function \code{decdiv} apportions trait-based diversity across the nodes of a phylogenetic tree using the algorithm defined in Pavoine et al. (2010). In that case the value at a given node is equal to the weight of a node times a measure of beta trait-based diversity between the lineages that descend from the node. The sum of all values attributed to the nodes of a phylogeny is equal to the total trait-based diversity of the species (tips of the phylogeny) as defined by Rao's quadratic entropy. 
In case of dichotomic trees, only two lineages descend from a given node. Here I consider a more general case where more than two lineages may descend from a node (polytomy). The beta trait-based diversity among the lineages that descend from a node is measured here as the average trait-based dissimilarity between any two of these descending lineages. With \code{option = 1}, the trait-based dissimilarity between two lineages is measured by Rao's DISC index (gamma diversity [average trait-based dissimilarity between any two species descending from the node] - alpha diversity [average trait-based dissimilarity between any two species descending from one of the lineages branched to the node]).   

In the present version of function \code{decdiv}, I have added other options. Options 2 and 3 code different ways of measuring trait-based differences between lineages, standardized between 0 and 1: with \code{option = 2}, the formula is (gamma - alpha) / (1 - alpha) * M / (M - 1), where gamma and alpha are defined above (for \code{option = 1}) and M is the number of lineages descending from the node; with \code{option = 3}, the formula is (gamma - alpha) / (1 - gamma) / (M - 1).

Options 4 and 5 decompose the result given by \code{option = 1}. \code{option = 4} returns gamma minus alpha (a measure of beta functional diversity between the lineages that descend from a node). \code{option = 5} returns the weights of the nodes (the summed relative abundance of its descending species or the relative number of descending species if presence/absence data are used). 
 
Values for \code{option} different from 1 needs that values in \code{dis} (the trait-based dissimilarities between species) are bounded between 0 and 1 if \code{formula = QE} or sqrt(2) if \code{formula = EDI}. If they are not bounded, the function \code{decdiv} will bound them using the maximum observed value in \code{dis}. 

The argument \code{ties.method} in \code{rtestdecdiv} allows you to explicitly take into account potential ties when ranking nodes according to their complexity or their distance to root (see Pavoine et al. 2010 for more details on the permutation test implemented in rtestdecdiv).
}
\references{
Champely, S. and Chessel, D. (2002) Measuring biological diversity using Euclidean metrics. \emph{Environmental and Ecological Statistics}, \bold{9}, 167--177.

Pavoine, S., Baguette, M., Bonsall, M.B. (2010) Decomposition of trait diversity among the nodes of a phylogenetic tree. \emph{Ecological Monographs}, \bold{80}, 485--507.

Rao, C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. \emph{Theoretical Population Biology}, \bold{21}, 24--43.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
\dontrun{
if(require(ape) && require(adephylo)){
data(ungulates)
ung.phy <- read.tree(text=ungulates$tre)
ung.phy$node.label
plot(ung.phy, show.node.label=TRUE)
# Regaring traits, we log-tranformed the first three traits 
# measuring volumes 
# and we standardized all the traits (mean=0; variance=1).
tab <- cbind.data.frame(afbw = log(ungulates$tab$afbw),
mnw = log(ungulates$tab$mnw), fnw = log(ungulates$tab$fnw),
ls = ungulates$tab$ls)
ung.tab0 <- data.frame(scalewt(tab))
ung.tab0 <- data.frame(scalewt(log(ungulates$tab)))
ung.pres <- rep(1, nrow(ung.tab0))
names(ung.pres) <- rownames(ung.tab0)
ung.dec1 <- decdiv(ung.phy, ung.pres, dist(ung.tab0), 
    option=1, formula = "EDI")
plot.decdiv(ung.dec1)

ung.dec2 <- decdiv(ung.phy, ung.pres, dist(ung.tab0), 
    option=2, formula = "EDI")
plot.decdiv(ung.dec2)

ung.dec3 <- decdiv(ung.phy, ung.pres, dist(ung.tab0), 
    option=3, formula = "EDI")
plot.decdiv(ung.dec3)

ung.dec4 <- decdiv(ung.phy, ung.pres, dist(ung.tab0), 
    option=4, formula = "EDI")
plot.decdiv(ung.dec4)

ung.dec5 <- decdiv(ung.phy, ung.pres, dist(ung.tab0), 
    option=5, formula = "EDI")
plot.decdiv(ung.dec5)
}
}
}
\keyword{models}
