\name{evodiss_family}
\alias{evodiss_family}
\alias{evodiss_ternaryplot}
\title{
A Family of Indices Dedicated to Pair-wise Phylogenetic Dissimilarities between Communities
}
\description{
The function \code{evodiss_family} was written thanks to function \code{\link{dist.binary}} of package ade4. Function \code{dist.binary} calculates specific compositional distances. The new function here replaces species with evolutionary units. It calculates Nipperess et al. (2010) parameters a, b, c, d (with incidence data), or A, B, C, D (with abundance data) and then use these parameters to compute pair-wise phylogenetic dissimilarities between communities.

The graphical function \code{evodiss_ternaryplot} displays Nipperess et al. (2010) parameters a, b, c (with incidence data), or A, B, C (with abundance data) on a ternary plot (see Koleff et al. 2003).
}
\usage{
evodiss_family(phyl, comm, method = NULL, abundance = TRUE,
squareroot = TRUE, diag = FALSE, upper = FALSE, tol = 1e-08)

evodiss_ternaryplot(phyl, comm, abundance = TRUE, 
tol = 1e-08, clabel = 1, ...)
}
\arguments{
  \item{phyl}{
an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase) or \code{hclust}.
}
  \item{comm}{
a data frame typically with communities as rows, species as columns and presence/absence or an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips.
}
  \item{method}{
a number between 1 and 10. If \code{NULL}, the choice is made with a console message. See details.
}
  \item{abundance}{
a logical indicating whether abundance data (if \code{TRUE}) or incidence data (presence/absence if \code{FALSE}) must be used.}
  \item{squareroot}{
a logical. First a similarity index (\emph{S}) is calculated (see details). Then if \code{squareroot=TRUE}, the dissimilarity is calculated as \eqn{\sqrt{1-S}}. If squareroot=FALSE, it is calculated as 1-\emph{S}. The square root is useful to obtain Euclidean properties for the PD-dissimilarity indices.}
  \item{diag}{
parameter (bolean) passed to function as.dist (R base).
}
  \item{upper}{
parameter (bolean) passed to function as.dist (R base).
}
  \item{tol}{
numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.
}
  \item{clabel}{
numeric: the size of the labels.}
   \item{\dots}{
other parameters can be added and passed to the function \code{triangle.plot} of the ade4 package (labels cannot be changed and are defined thanks to the row names of \code{comm}).
}
}
\details{
The function was written thanks to function \code{dist.binary} of package ade4. Function \code{dist.binary} calculates specific compositional distances. The new function here replaces species with evolutionary units. It calculates Nipperess et al. (2010) parameters a, b, c, d (with incidence data), or A, B, C, D (with abundance data). Then, the parameters are combined thanks to one out of ten methods.

\code{1}: Jaccard index (1901) S3 coefficient of Gower and Legendre (1986) s1 = a / (a+b+c).

\code{2}: Simple matching coefficient of Sokal and Michener (1958) S4 coefficient of Gower and Legendre (1986) s2 = (a+d) / (a+b+c+d).

\code{3}: Sokal and Sneath(1963) S5 coefficient of Gower and Legendre (1986) s3 = a / (a + 2(b + c)).

\code{4}: Rogers and Tanimoto (1960) S6 coefficient of Gower and Legendre (1986) s4 = (a + d) / (a + 2(b + c) +d).

\code{5}: Dice (1945) or Sorensen (1948) S7 coefficient of Gower and Legendre (1986) s5 = 2a / (2a + b + c).

\code{6}: Hamann coefficient S9 index of Gower and Legendre (1986) s6 = (a - (b + c) + d) / (a + b + c + d).

\code{7}: Ochiai (1957) S12 coefficient of Gower and Legendre (1986) s7 = a / sqrt((a + b)(a + c)).

\code{8}: Sokal and Sneath (1963) S13 coefficient of Gower and Legendre (1986) s8 = ad / sqrt((a + b)(a + c)(d + b)(d + c)).

\code{9}: Phi of Pearson S14 coefficient of Gower and Legendre (1986) s9 = (ad - bc) / sqrt((a + b)(a + c)(d + b)(d + c)).

\code{10}: S2 coefficient of Gower and Legendre (1986) s10 = a / (a + b + c + d).
}
\value{
Function \code{evodiss_family} returns an object of class \code{dist} containing the PD-dissimilarities (phylogenetic dissimilarities) between communities.

Function \code{evodiss_ternaryplot} returns a graph.
}
\references{
The methodologies are presented in 
Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 

They gather in a common framework and extend earlier work introduced in  
Koleff, P., Gaston, K.J., Lennon, J.J. (2003) Measuring beta diversity for presence-absence data. \emph{Journal of Animal Ecology}, \bold{72}, 367--382.

Nipperess, D.A., Faith, D.P., Barton, K.(2010) Resemblance in phylogenetic diversity among ecological assemblages. \emph{Journal of Vegetation Science}, \bold{21}, 809--820.

Dice, L.R. (1945) Measures of the amount of ecologic association between species. \emph{Ecology} \bold{26}, 297--302.

Gower, J.C., Legendre, P. (1986) Metric and Euclidean properties of dissimilarity coefficients. \emph{Journal of Classification}, \bold{3}, 5--48.

Jaccard, P. (1901) Etude comparative de la distribution florale dans une portion des Alpes et des Jura. \emph{Bulletin de la Societe Vaudoise des Sciences Naturelles}, \bold{37}, 547--579.

Ochiai, A. (1957) Zoogeographic studies on the soleoid fishes found in Japan and its neighbouring regions. \emph{Bulletin of the Japanese Society of Scientific Fisheries}, \bold{22}, 526--530.

Rogers, J.S. and Tanimoto, T.T. (1960) A computer program for classifying plants. \emph{Science}, \bold{132}, 1115--1118.

Sokal, R.R. and Michener, C.D. (1958) A Statistical Method for Evaluating Systematic Relationships. \emph{The University of Kansas Science Bulletin}, \bold{38}, 1409--1438.

Sokal, R.R. and Sneath, P.H.A. (1963) \emph{Principles of numerical taxonomy}. San Francisco: W. H. Freeman.

Sorensen, T. 1948. A method of establishing groups of equal amplitude in plant sociology based on similarity of species content. \emph{Kongelige Danske Videnskabernes Selskabs Biologiske Skrifter}, \bold{5}, 1--34.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evodiss}}
}
\examples{
data(batcomm)
phy <- read.tree(text=batcomm$tre)
ab <- batcomm$ab[,phy$tip.label]

# PD-dissimilarity indices that use Nipperess et al. (2010)
# parameters can be obtained thanks to function evodiss_family. 
# For example, with incidence data, 
# indices evoDJaccard, evoDSorensen, and evoDOchiai 
# (supplementary Appendix 1 in Pavoine 2016) 
# can be obtained as follows:

evodiss_family(phy, ab, method=1, abundance=FALSE) # Jaccard
evodiss_family(phy, ab, method=5, abundance=FALSE) # Sorensen
evodiss_family(phy, ab, method=7, abundance=FALSE) # Ochiai

# With abundance data, indices evoDTJ, evoDTS, evoDTO 
# (supplementary Appendix 1 in Pavoine 2016) 
# can be obtained as follows:

evodiss_family(phy, ab, method=1) # evoDTJ
evodiss_family(phy, ab, method=5) # evoDTS
evodiss_family(phy, ab, method=7) # evoDTO

# Ternary plots can be obtained for Nipperess et al. (2010) 
# parameters a, b, c (incidence data) 
# (see Supplementary material Appendix 4 in Pavoine 2016):

evodiss_ternaryplot(phy, ab, abundance = FALSE)

# and for Nipperess et al. (2010) parameters A, B, C 
# (abundance data):

evodiss_ternaryplot(phy, ab, abundance = TRUE)

# The ternary plots can be adjusted thanks 
# to parameters of function triangle.plot (package ade4). 
# For example, full triangles can be obtained as follows 
# (previous graphs were zoomed on the smallest principal 
# equilateral triangle that contained the points, 
# as indicated by the embedded close grey triangle 
# at the left-hand corner of ternary plot given above):

evodiss_ternaryplot(phy, ab, abundance = FALSE, scale=FALSE, show.position=FALSE)
# Incidence data

evodiss_ternaryplot(phy, ab, abundance = TRUE, scale=FALSE, show.position=FALSE)
# abundance data
}
\keyword{models}

