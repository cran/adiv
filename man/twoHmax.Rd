\name{twoHmax}
\alias{twoHmax}
\title{
Maximization of the Diversity Index TwoH
}
\description{
The R function \code{twoHmax} maximizes function 2H (Pavoine and Izsak 2014) for a given matrix \bold{C} of (functional or phylogenetic) similarities between species. It is based on function \code{divcmax} of the ade4 package in R. As function \code{divcmax} of package ade4, function \code{twoHmax} uses an optimization technique based on Rosen's projection gradient algorithm and is verified using the Kuhn-Tucker conditions.
}
\usage{
twoHmax(C, epsilon = 1e-08, smooth = TRUE, comment = FALSE)
}
\arguments{
  \item{C}{
a matrix that contains measures of the chosen intraspecific components (see functions \code{CFprop}, \code{multiCFprop}, \code{CFbinary}, \code{multiCFbinary}) on the diagonal and measures of interspecific components off diagonal. These interspecific components reflect functional or phylogenetic similarities among species.
}
  \item{epsilon}{
a numeric tolerance threshold: a frequency is non null if it is higher than epsilon.
}
  \item{smooth}{
a logical value: if \code{TRUE}, the estimated positive proportions are checked and re-estimated with exact matrix products.}
  \item{comment}{
a logical value indicating whether or not comments on the optimization technique should be printed.}
}
\value{
A list of two objects:
\item{value}{the maximum value of index 2H (see function \code{qHdiv});}
\item{vector}{a data frame with the vector pmax that maximizes index 2H (see function \code{qHdiv}).}
}
\references{
Pavoine, S. and Izsak, J. (2014) New biodiversity measure that includes consistent interspecific and intraspecific components. \emph{Methods in Ecology and Evolution}, \bold{5}, 165--172.}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}

The code is a modification of function divcmax of package ade4 written by Stephane Champely.
}
\seealso{
\code{\link{qHdiv}}, \code{\link{CFprop}}, \code{\link{CFbinary}}
}
\examples{
\dontrun{
if(require(ape)){

tre <-"((((((sA:4,sB:1):1,sC:3):2,((sD:2,sE:1):1,sF:1):2):1,sG:7):1,sH:1):3,(sI:2,sJ:1):2):0;"
#The number of tips is kept in parameter n:
n<-10
# Next we need to obtain matrix CP. 
phyape <- read.tree(text = tre)
plot(phyape)
CP <- vcv(phyape)
WP <- diag(diag(CP))
# With this particular illustration, a maximizing vector, 
# for 2H used with CP, that does not contain any zero can be found. 
# This maximizing vector can thus be obtained directly, 
# instead of being estimated. Two equivalent equations 
# have been given to obtain the maximizing vector in 
# Appendix S1 of Pavoine and Izsak (2014). 
# We use the first one below
Pmax<-(solve(CP^2)\%*\%diag(CP))/sum(solve(CP^2)\%*\%diag(CP))
Pmax

# The second equation equivalently provides
Z <- ((diag(1/sqrt(diag(CP))))\%*\%CP\%*\%(diag(1/sqrt(diag(CP)))))^2
Pmax<-(solve(WP)\%*\%solve(Z)\%*\%rep(1,n))/sum(solve(WP)\%*\%solve(Z))
Pmax

# Applied to our case study, the function twoHmax provides good approximations
twoHmax(CP)

# Redundancy among variables:
data(rhone, package="ade4")
V <- rhone$tab
# First consider the covariances among the variables:
C <- cov(V)
# A vector that maximizes 2H applied to C is estimated as follows:
pmax_covariances <- twoHmax(C)$vector
dotchart(as.matrix(pmax_covariances))

# If we apply 2H only to the diagonal matrix with the variances 
# of the variables, the vector that maximizes 2H is:
W <- diag(diag(C))
rownames(W)<-colnames(W)<-rownames(C)
pmax_variances <- twoHmax(W)$vector
dotchart(as.matrix(pmax_variances))

# If C contains the correlations among variables, 
# a vector that maximizes 2H applied to C is estimated as follows:
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
# the effective number of variables is lower.
C <- cor(V)
equalproportions <- cbind.data.frame(rep(1/ncol(C), ncol(C)))
names(equalproportions) <- "equalprop"
equalproportions <- t(equalproportions)
qHdiv(equalproportions, C)

# When considering the covariances among species, 
# instead of the correlations, the effective number 
# of variables is even lower, indicating also an imbalance 
# in the variances of the variables.
C <- cov(V)
qHdiv(equalproportions, C)

}
}
}
\keyword{models}
