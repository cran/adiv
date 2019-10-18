\name{qHdiv}
\alias{qHdiv}
\title{
Biodiversity Measure that Includes Consistent Interspecific and Intraspecific Components
}
\description{
Function \code{qHdiv} calculates the parametric diversity index developed by Pavoine and Izsak (2014)
}
\usage{
qHdiv(comm, C, q = 2)
}
\arguments{
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and an index of abundance as entries.
}
  \item{C}{
a matrix that contains measures of the chosen intraspecific components (as defined in Pavoine and Izsak 2014) on the diagonal and measures of interspecific components off diagonal. These interspecific components reflect functional or phylogenetic similarities among species.
}
  \item{q}{
a positive numeric for the value of the parameter \code{q}.
}
}
\value{
a vector with the diversity in each of the communities (same order as in \code{comm}).
}
\references{
Pavoine, S., Izsak, J. (2014) New biodiversity measure that includes consistent interspecific and intraspecific components. \emph{Methods in Ecology and Evolution}, \bold{5}, 165--172.}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{twoHmax}}, \code{\link{CFprop}}, \code{\link{CFbinary}}
}
\examples{
\dontrun{
if(require(ape)){

# First case study on phylogenetic diversity:
# Below is a phylogenetic tree in newick format
tre <-"((((((sA:4,sB:1):1,sC:3):2,((sD:2,sE:1):1,sF:1):2):1,sG:7):1,sH:1):3,(sI:2,sJ:1):2):0;"
#The number of tips is kept in parameter n:
n<-10
# Next we need to obtain matrix CP (see Pavoine and Izsak 2014). 
phyape <- read.tree(text = tre)
plot(phyape)
CP <- vcv(phyape)
WP <- diag(diag(CP))
# With this particular illustration, a maximizing vector, 
# for 2H applied to CP, that does not contain any zero can be found. 
# This maximizing vector can thus be obtained directly, 
# instead of being estimated. 
# Two equivalent equations have been given to obtain 
# the maximizing vector in Appendix S1 of 
# Pavoine and Izsak (2014). 
# We use the first one below:

Pmax<-(solve(CP^2)\%*\%diag(CP))/sum(solve(CP^2)\%*\%diag(CP))
Pmax

# The second equation equivalently provides:

Z <- ((diag(1/sqrt(diag(CP))))\%*\%CP\%*\%(diag(1/sqrt(diag(CP)))))^2
Pmax<-(solve(WP)\%*\%solve(Z)\%*\%rep(1,n))/sum(solve(WP)\%*\%solve(Z))
Pmax

# Applied to our case study, the function twoHmax 
# provides a good approximation of the maximizing vector:

twoHmax(CP)

# Second case study on the redundancy among variables:
data(rhone, package="ade4")
V <- rhone$tab
# First consider the covariances among the variables:
C <- cov(V)
# A vector that maximizes 2H applied to C is estimated 
# as follows:
pmax_covariances <- twoHmax(C)$vector
dotchart(as.matrix(pmax_covariances))

# If we apply 2H only to the diagonal matrix with the variances 
# of the variables, the vector that maximizes 2H is:
W <- diag(diag(C))
rownames(W)<-colnames(W)<-rownames(C)
pmax_variances <- twoHmax(W)$vector
dotchart(as.matrix(pmax_variances))

# If C contains the correlations among variables, 
# a vector that maximizes 2H applied to C is estimated 
# as follows:
C <- cor(V)
pmax_correlations <- twoHmax(C)$vector
dotchart(as.matrix(pmax_correlations))

# By attributing equal weights to the variables, 
# 2H applied to the correlation matrix measures 
# the number of effective variables:
# from 0 if all variables are completely correlated 
# with each other to n if they are not correlated. 
# Similarly, by attributing equal weights to the variables, 
# 2H applied to the covariance matrix measures 
# the effective amount of variation:
# from 0 if all variables are completely correlated 
# with each other to n if they are not correlated 
# and have similar variances.

#Even if the data set contains 15 variables, 
# the effective number of variables is lower:
C <- cor(V)
equalproportions <- cbind.data.frame(rep(1/ncol(C), ncol(C)))
names(equalproportions) <- "equalprop"
equalproportions <- t(equalproportions)
qHdiv(equalproportions, C)

# When considering the covariances, instead of the correlations, 
# the effective number of variables is even lower,
# indicating also an imbalance in the variances 
# of the variables.
C <- cov(V)
qHdiv(equalproportions, C)

}
}
}
\keyword{models}
