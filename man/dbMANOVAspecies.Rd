\name{dbMANOVAspecies}
\alias{dbMANOVAspecies}
\alias{dbMANOVAspecies_pairwise}
\alias{summary.dbMANOVAspecies_pairwise}
\title{Indicator species analysis in the framework of multivariate analysis of variance}
\description{
"Consider N plots distributed among K groups; plots are composed of species whose abundances within each plot are known. The QT, QB and QW statistics defined in [Ricotta et al. 2021] aim at evaluating the average difference, in terms of species identity and abundance, between any two plots (QT), between two plots within a group (QW), and the gap between QT and QW (=QB) due to compositional differences between groups of plots. The function \code{dbMANOVAspecies} calculates QT, QB and QW as in [Ricotta et al. 2021] and the species-centered components QTs, QBs and QWs. It also calculates the SES values (equation 7 [in Ricotta et al. 2021]) associated with QB and QBs and allows tests for the significance of the SES values (H0 = the SES value is similar as that expected by randomly permuting plots among groups of plots and H1 = the SES value of is greater than that expected by randomly permuting plots among groups of plots). The function \code{dbMANOVAspecies_pairwise} complements function \code{dbMANOVAspecies} by performing post-hoc tests for all pairs of groups. It must be executed in the same R environment (workspace) as function \code{dbMANOVAspecies}. A third function named \code{summary.dbMANOVAspecies_pairwise} provides a short summary of the results of function \code{dbMANOVAspecies_pairwise} (with SES and P values)" (Ricotta et al. 2021, Appendix 3).
}
\usage{
dbMANOVAspecies(comm, groups, nrep = 999, method = c("Euclidean", "Manhattan", 
    "Canberra", "BrayCurtis"), global = TRUE, species = TRUE, padjust = "none", 
    tol = 1e-8)
    
    dbMANOVAspecies_pairwise(dbobj, signif = TRUE, salpha = 0.05, nrep = NULL)
    
    \method{summary}{dbMANOVAspecies_pairwise}(object, DIGITS = 3, \dots)
}
\arguments{
\item{comm}{a matrix of N plots \eqn{\times} S species containing the relative or absolute abundance of all species. Columns are species and plots are rows.}
\item{groups}{a vector of characters or a factor with the names of the groups associated with plots. Names of groups must be listed in the same order as plots in object comm. For example, the first value of groups gives the name of the group for the first plot (first row in object comm).}
\item{nrep}{a numeric that gives the number of permutations to be done. nrep can be set to NULL in function \code{dbMANOVAspecies_pairwise}, in which case the number of permutations used to built object dbobj (first argument of function \code{dbMANOVAspecies_pairwise}) is used by default.}
\item{method}{a string, one of \code{"Euclidean"} (to perform the Euclidean distance), \code{"Manhattan"} (to select the Manhattan distance), \code{"Canberra"} (to perform the Canberra distance), \code{"BrayCurtis"} (to opt for the BrayCurtis distance).}
\item{global}{a logical. If \code{TRUE} the global test (using the SES of QB) is performed. The global test can be performed in addition to one test per species.}
\item{species}{a logical. If \code{TRUE} one test per species is performed. If \code{species = FALSE}, then argument global is ignored and the global test is performed.}
\item{padjust}{a method of correction for multiple tests one of p.adjust.methods (see \code{?p.adjust.methods} for possible choices); choose \code{"none"} if you don't need any correction for multiple tests. The argument \code{padjust} is ignored if \code{species = FALSE}.}
\item{tol}{a numeric tolerance threshold. Any value between \code{-tol} and \code{tol} will be considered as equal to zero.}
\item{dbobj}{an object of class \code{"dbMANOVAspecies"} obtained with function \code{dbMANOVAspecies}.}
\item{signif}{a logical. If \code{TRUE}, the post-hoc tests are performed only for species that were associated with significant tests in object dbobj.}
\item{salpha}{a numeric. The level of significance (nominal alpha error) for P values (must be between 0 and 1). Ignored if \code{signif = FALSE}.}
\item{object}{an object of class \code{"dbMANOVAspecies_pairwise"} obtained with function \code{dbMANOVAspecies_pairwise}.}
\item{DIGITS}{integer indicating the number of decimal places to retain when displaying the results. If NULL, all decimals are retained.}
\item{\dots}{further arguments passed to or from other methods.}
}
\value{
\code{dbMANOVAspecies} returns a list with the following objects:

-\code{observations}: a vector or a data frame with the observed values of statistics (QT, QB and QW) (if global = TRUE) and/or the contribution each species has to these statistics (if species = TRUE);

-\code{test}: an object of class randtest or krandtest with the results of the test for the global differences, in terms of species composition, between groups of plots (if global = TRUE) and/or for the contribution each species has in these differences (if species = TRUE). \code{dbMANOVAspecies_pairwise} provides in a list the same objects as function \code{dbMANOVAspecies} but for each pair of groups (see examples below).

If only a global test was performed, \code{summary.dbMANOVAspecies_pairwise} provides a data frame with two rows named "SES" for the SES values (of the QB statistic) and "pvalue" for the P values and as many columns as there are pairs of groups of plots. Columns are named according to the two groups that are compared. For example, if there are two groups named Group1 and Group2, they are compared in column named "Group1:Group2". Else, it provides a data frame with species names as rows (in addition to a row named GLOBAL if a global test was also performed) and as columns the SES, P values, and if relevant adjusted P values, for each combination of two groups. Consider that Group1 and Group2 are the names of two groups. The column names that correspond to the comparison between Group1 and Group2 are written Group1:Group2.SES (for the SES values), Group1:Group2.pvalue (for the P values), and Group1:Group2.adj.pvalue (for the P values adjusted after correction for multiple tests, if a correction was done). NAs (= missing values) may be present in the data frame for species that are absent from two compared groups of plots. NAs may also be present in the data frame if signif = TRUE was used in function \code{dbMANOVAspecies_pairwise} for species that were associated to a non-significant test when performing function \code{dbMANOVAspecies}.

}
\references{
Ricotta, C., Pavoine, S., Cerabolini, B.E.L., Pillar, V.D. (2021) A new method for indicator species analysis in the framework of multivariate analysis of variance. \emph{Journal of Vegetation Science}, \bold{32}, e13013. \doi{10.1111/jvs.13013}
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
\dontrun{
data(RutorGlacier)
Qspecies <- dbMANOVAspecies(RutorGlacier$Abund, RutorGlacier$Fac, 
    nrep=9999, global=TRUE)
Qspecies_adj <- dbMANOVAspecies(RutorGlacier$Abund, RutorGlacier$Fac, 
    nrep=9999, global=FALSE, padj = "BY")
# In Qspecies and Qspecies_adj, Column "Std.Obs" contains the SES values.
# Now for all species that showed significant compositional difference among 
# the three groups of plots (with a nominal alpha error of 0.05), we can test for 
# pairwise differences among all pairs of groups thanks to function 
# dbMANOVAspecies_pairwise as shown below:

# without any correction for multiple tests
Qspeciespairwise <- dbMANOVAspecies_pairwise(Qspecies)
summary(Qspeciespairwise)
# NAs are present in the data frame above for species that were associated 
# to non-significant test in object Qspecies (where tests were done over all groups);
# and also, for species that are absent from the two compared groups of plots 
# (e.g. species Adenostyles leucophylla, in mid- and late-successional stages).

# with correction for multiple test
Qspeciespairwise_adj <- dbMANOVAspecies_pairwise(Qspecies_adj)
summary(Qspeciespairwise_adj)
# Here again, NAs are present in the data frame above for species that were 
# associated to non-significant test in object Qspecies_adj 
# (where tests were done over all groups); and also, for species that are absent 
# from the two compared groups of plots (e.g. species Adenostyles leucophylla, 
# in mid- and late-successional stages).

# See Ricotta et al. 2021 Appendix 3 for details.
}
}
\keyword{models}

