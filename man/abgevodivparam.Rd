\name{abgevodivparam}
\alias{abgevodivparam}
\alias{plot.abgevodivparam}
\title{
Apportionment of Parametric Indices of Phylogenetic Diversity
}
\description{
Function \code{abgevodivparam} calculates alpha, beta and gamma components of phylogenetic diversity using parametric indices derived from Tsallis (HCDT) and Hill compositional indices. Alpha is for within-community diversity, beta for between-community diversity and gamma for the diversity of all combined communities.
}
\usage{
abgevodivparam(phyl, comm, w = c("evoab", "even", "speciesab"),
method = c("hillCJC", "hillR", "tsallis"), q = 2, 
option = c("multiplicative", "additive", "proportional", 
"C", "U", "V", "S", "Renyi"), tol = 1e-08)

\method{plot}{abgevodivparam}(x, legend = TRUE, 
legendposi = "topright", type = "b", 
col = if (is.numeric(x)) NULL else 1:nrow(x$div), 
lty = if (is.numeric(x)) NULL else rep(1, nrow(x$div)), 
pch = if (is.numeric(x)) NULL else 1:nrow(x$div), 
ylim1 = range(x$div[c("Alpha", "Gamma"), ]), ylim2 = NULL, ...)
}
\arguments{
  \item{phyl}{an object inheriting the class \code{phylo} (see package ape), \code{phylo4} (see package phylobase), or \code{hclust}.}
  \item{comm}{a data frame or a matrix typically with communities as rows, species as columns and an index of abundance as entries. Species should be labeled as in the phylogenetic tree where they are the tips.}
  \item{w}{either a numeric vector giving weights for communities (same order as in comm), or a code: one of \code{"even"}, \code{"evoab"}, and \code{"speciesab"}. If several codes are given, only the first one is used. See details.}
  \item{method}{a string with one of the following codes: \code{"tsallis"}, \code{"hillR"}, or \code{"hillCJC"}. See details.}
  \item{q}{a vector with nonnegative value(s) for parameter \code{q}. See details.}
  \item{option}{a string code: either \code{"multiplicative"}, \code{"additive"} or \code{"proportional"}. If several codes are given, only the first one is used. Only, with \code{method}=\code{"hillCJC"}, other options are possible: \code{"C"}, \code{"U"}, \code{"V"}, \code{"S"}, \code{"Renyi"}. See details.}
  \item{tol}{numeric tolerance threshold: values between -\code{tol} and \code{tol} are considered equal to zero.}
  \item{x}{an object of class \code{abgevodivparam} obtained with function \code{abgevodivparam}.}
  \item{legend}{a logical. If TRUE a legend is given with the colour, the type of line (etc.) used to define the diversity curve of each diversity level (gamma, alpha, beta).}
  \item{legendposi}{a string that gives the position of the legend to be passed to function \code{legend} of the base of R.}
  \item{type}{a string to be passed to the graphic argument \code{type} of functions \code{plot} and \code{lines} used to draw the diversity curve of each diversity level (gamma, alpha, beta).}
  \item{col}{vector of colours to be passed to the graphic argument \code{col} of functions \code{plot} and \code{lines} to define the colour of the diversity curve of each diversity level (gamma, alpha, beta, in that order).}
  \item{lty}{vector of types of line (plain, broken etc.) to be passed to the graphic argument \code{lty} of functions \code{plot} and \code{lines} used to draw the diversity curve of each diversity level (gamma, alpha, beta, in that order).}
  \item{pch}{vector of types of point (open circle, close circle, square etc.) to be passed to the graphic argument \code{pch} of functions \code{plot} and \code{lines} used to draw the diversity level (gamma, alpha, beta, in that order).}
  \item{ylim1}{a vector with two numerics indicating the range to be used to display alpha and gamma diversity.}
  \item{ylim2}{a vector with two numerics indicating the range to be used to display beta diversity.}
  \item{\dots}{other arguments can be added and passed to the functions \code{plot} and \code{lines} used to draw the graphic.}
}
\details{
Consider a phylogenetic tree \emph{T}, \eqn{b_T}{b_T} the set of branches in \emph{T}, \emph{k} a branch, \eqn{L_k}{L_k} the length of branch \emph{k}, \emph{j} a community (\emph{j}=1,...,\emph{m}), \eqn{a_{jk}}{a_jk} the abundance associated with branch \emph{k} in community \emph{j} (sum of abundance of all species descending from the branch). \emph{q} is the parameter that increases with the importance given to abundant species compared to rare species in diversity.

The methods available are:
\code{tsallis} (decomposition of Tsallis or HCDT entropy (Harvda and Charvat 1967; Daroczy 1970; Tsallis 1988) into alpha, beta, gamma components adapted here to phylogenetic diversity):

\deqn{^q\gamma_{evoTsallis}=\left[1-\sum_{k \in b_T} L_k \left(\sum_{j=1}^m w_j \frac{a_{jk}}{\sum_{k \in b_T} L_k a_{jk}}\right)^q \right]/(q-1)}{q\gamma_evoTsallis=[1-sum_k L_k (sum_j w_j a_jk/(sum_k L_k a_jk))^q]/(q-1)}

\deqn{^q\alpha_{evoTsallis}=\sum_{j=1}^m w_j \left[1-\sum_{k \in b_T} L_k \left( \frac{a_{jk}}{\sum_{k \in b_T} L_k a_{jk}}\right)^q\right]/(q-1)}{q\alpha_evoTsallis=sum_j w_j [1-sum_k L_k (a_jk/(sum_k L_k a_jk))^q]/(q-1)}

\code{hillR} (Routledge decomposition of Hill diversity into alpha, beta, gamma components adapted hete to phylogenetic diversity):

\deqn{^q\gamma_{evoHill}=\left[\sum_{k \in b_T} L_k \left(\sum_{j=1}^m w_j \frac{a_{jk}}{\sum_{k \in b_T} L_k a_{jk}}\right)^q \right]^{1/(1-q)}}{q\gamma_evoHill=[sum_k L_k (sum_j w_j a_jk/(sum_k L_k a_jk))^q]^(1/(1-q))}

\deqn{^q\alpha_{evoHill-R}=\left[\sum_{j=1}^m w_j \sum_{k \in b_T} L_k \left( \frac{a_{jk}}{\sum_{k \in b_T} L_k a_{jk}}\right)^q\right]^{1/(1-q)}}{q\alpha_evoHill-R=[sum_j w_j sum_k L_k (a_jk/(sum_k L_k a_jk))^q]^(1/(1-q))}

\code{hillCJC} (Chiu et al. (2014) decomposition of phylogenetic diversity into alpha, beta, gamma components, see Supplementary material Appendix 2 in Pavoine (2016) for a justification of the formulas):

\deqn{^q\gamma_{evoHill}=\left[\sum_{k \in b_T} L_k \left(\sum_{j=1}^m w_j \frac{a_{jk}}{\sum_{k \in b_T} L_k a_{jk}}\right)^q \right]^{1/(1-q)}}{q\gamma_evoHill=[sum_k L_k (sum_j w_j a_jk/(sum_k L_k a_jk))^q]^(1/(1-q))}

\deqn{^q\alpha_{evoHill-CJC}=\frac{1}{m}\left[\sum_{k \in b_T} L_k \sum_{j=1}^m (w_j)^q \left( \frac{a_{jk}}{\sum_{k \in b_T} L_k a_{jk}}\right)^q\right]^{1/(1-q)}}{q\alpha_evoHill-CJC=(1/m)*[sum_k L_k sum_j (w_j)^q (a_jk/(sum_k L_k a_jk))^q]^(1/(1-q))}

Then option \code{"additive"} calculates \eqn{\beta} diversity as \eqn{\gamma-\alpha}.
Option \code{"proportional"} calculates \eqn{\beta} as \eqn{(\gamma-\alpha)/\gamma}.
Option \code{"multiplicative"} calculates \eqn{\beta} diversity as \eqn{\gamma/\alpha}.
Only for \code{method}=\code{"hillCJC"}, options \code{"C"}, \code{"U"}, \code{"V"}, \code{"S"}, use the multiplicative option and also calculate one of the transformations introduced by Chiu et al. (2014): indices \eqn{1-C_{qm}}{1-C_qm}, \eqn{1-U_{qm}}{1-U_qm}, \eqn{1-V_{qm}}{1-V_qm}, and \eqn{1-S_{qm}}{1-S_qm}, respectively. \code{"Renyi"} is the \eqn{^qevoD_{Renyi}}{qevoD_Renyi} index introduced in Pavoine (2016), see also Supplementary material Appendix 1 in Pavoine (2016).

The weights of the sites (argument \code{w}) can be \code{"even"} (even weights), \code{"evoab"} (proportional to the summed abundances of all evolutionary units), or \code{"speciesab"} (proportional to the summed abundances of all species). Note that if the phylogenetic tree is ultrametric (the distance from any species to the root is constant), then options \code{"evoab"} and \code{"speciesab"} are equivalent.
}
\value{
If only one value of \code{q} is given, abgevodivparam returns a vector with alpha, beta, and gamma diversities.
If more than one value of \code{q} is given, it returns a list of two objects: 
\item{q}{the numeric vector of values for \code{q}.}
\item{div}{a data frame with alpha, beta, gamma calculated for all values of \code{q}.}
Only if \code{method}=\code{"hillCJC"} and \code{option}= \code{"C"}, \code{"U"}, \code{"V"}, \code{"S"}, or \code{"Renyi"}, the index \eqn{1-C_{qm}}{1-C_qm} (for \code{"C"}), \eqn{1-U_{qm}}{1-U_qm} (for \code{"U"}), \eqn{1-V_{qm}}{1-V_qm} (for \code{"V"}), \eqn{1-S_{qm}}{1-S_qm} (for \code{"S"}) or the Renyi transformation (see above, for \code{"Renyi"} is also provided in the \code{div} data frame under the name "transformed.beta".

The function \code{plot.abgevodivparam} returns a graphic.
}
\references{
The methodologies and scripts were presented in 

Pavoine, S. (2016) A guide through a family of phylogenetic dissimilarity measures among sites. \emph{Oikos}, \bold{125}, 1719--1732. 

using earlier work by

Daroczy, Z. (1970) Generalized information functions. \emph{Information and Control}, \bold{16}, 36--51.

Havrda, M., Charvat F. (1967) Quantification method of classification processes: concept of structural alpha-
entropy. \emph{Kybernetik}, \bold{3}, 30--35

Hill, M.O. (1973) Diversity and evenness: a unifying notation and its consequences. \emph{Ecology}, \bold{54}, 427--432.

Routledge, R.D. (1979) Diversity indices: which ones are admissible? \emph{Journal of Theoretical Biology}, \bold{76}, 503--515.

Rao, C.R. (1986) Rao's axiomatization of diversity measures. In: Kotz S, Johnson NL, editors. \emph{Encyclopedia of Statistical Sciences}. New York: Wiley and Sons. pp. 614--617.

Chiu, C.-H., Jost, L., Chao, A. (2014) Phylogenetic beta diversity, similarity, and differentiation measures based on Hill numbers. \emph{Ecological Monographs}, \bold{84}, 21--44.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\seealso{
\code{\link{evodiss}}, \code{\link{divparam}}, \code{\link{evodivparam}}
}
\examples{
\dontrun{
if(require(ape)){
data(batcomm)
phy <- read.tree(text=batcomm$tre)
ab <- batcomm$ab[,phy$tip.label]
abgevodivparam(phy, ab)
plot(abgevodivparam(phy, ab))
abgevodivparam(phy, ab, q=0:4)
plot(abgevodivparam(phy, ab, q=0:4))
}
}
}
\keyword{models}
