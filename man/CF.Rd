\name{CF}
\alias{CFbinary}
\alias{multiCFbinary}
\alias{CFprop}
\alias{multiCFprop}
\title{
Matrices of Intra- and Inter-Specific Similarities
}
\description{
The function \code{CFprop} calculates the matrix CF of intra- (diagonal) and inter-specific (off-diagonal) similarities as defined in the main text of Pavoine and Izsak (2014), and matrix CwF as defined in Appendix S3 of Pavoine and Izsak (2014) for weighting functional attributes. The function \code{multiCFprop} calculates matrices CwmF1, CwmF2, CwmF3 when several functional traits are considered (Appendix S3 of Pavoine and Izsak 2014). Traits and the attributes of the traits can be weighted. These two functions consider functional traits expressed as proportion (compositional) vectors. The functions \code{CFbinary} and \code{multiCFbinary} are the equivalents of \code{CFprop} and \code{multiCFprop} when traits are expressed as binary vectors as shown in Appendix S3 of Pavoine and Izsak (2014).}
\usage{
CFbinary(df, wA = rep(1, ncol(df)))

multiCFbinary(Ktab, w.attributes = lapply(Ktab, function(x) rep(1, ncol(x))), 
w.traits = rep(1/length(Ktab), length(Ktab)), 
labels = rownames(Ktab[[1]]), solution = c(2, 1))

CFprop(df, wA = rep(1, ncol(df)))

multiCFprop(Ktab, w.attributes = lapply(Ktab, function(x) rep(1, ncol(x))), 
w.traits = rep(1/length(Ktab), length(Ktab)), 
labels = rownames(Ktab[[1]]), solution = c(2, 1))
}
\arguments{
  \item{df}{a data frame or a matrix with species (or any entities of interest) as rows, functional attributes as columns, and proportions (for \code{CFprop}) or 0/1 values (for \code{CFbinary}) as entries. The sum of the columns will be standardized to equal 1 in function \code{CFprop} but not in \code{CFbinary}.}
  \item{wA}{a vector of weights that should be given to the attributes (same order as the columns of df).}
  \item{Ktab}{a list of data frames, each of which represents a trait. For a given trait, the data frame should have species (or any entities of interest) as rows, functional attributes as columns, and proportions (for \code{multiCFprop}) or 0/1 values (for \code{multiCFbinary}) as entries. The sum of the columns will be standardized to equal 1 in function \code{multiCFprop} but not in \code{multiCFbinary}.}
  \item{w.attributes}{a list of weights that should be given to the attributes of each trait. Traits should be in the same order as they appear in the list of tables \code{Ktab}. The attributes of a trait should be ordered as the columns of the corresponding data frame in \code{Ktab}.}
  \item{w.traits}{a numeric vector of weights that should be given to the traits (same order as the tables of \code{Ktab}).}
  \item{labels}{a vector of strings that gives the names of the species (or the other entities of interest; same order as the rows of all tables in \code{Ktab}).}
  \item{solution}{a numeric value (either 1 or 2) that indicates which equations are used to summarize the information given by several traits among the 2 approaches given in Appendix S3 of Pavoine and Izsak (2014) page 9. If a vector is given, only the first value of the vector is considered.}
}
\value{
A matrix with nonnegative values
}
\references{
Pavoine, S. and Izsak, J. (2014) New biodiversity measure that includes consistent interspecific and intraspecific components. \emph{Methods in Ecology and Evolution}, \bold{5}, 165--172.}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{qHdiv}}
}
\keyword{models}
