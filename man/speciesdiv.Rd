\name{speciesdiv}
\alias{speciesdiv}
\title{
Indices of Species Diversity
}
\description{
The function \code{speciesdiv} calculates diversity indices that rely on relative or absolute species abundance. 
}
\usage{
speciesdiv(comm, method = "full", tol = 1e-8)
}
\arguments{
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and presence/absence or abundance as entry. Note that with presence/absence (0/1) data, only species richness will be calculated correctly.
}
  \item{method}{
a string or a vector of strings: one or several of "richness", "GiniSimpson", "Simpson", "Shannon", "Margalef", "Menhinick", "McIntosh", "full". See details.
}
  \item{tol}{
a tolerance threshold (a value between -\code{tol} and \code{tol} is considered equal to zero)}
}
\value{
Function \code{speciesdiv} returns a matrix with communities as rows and the diversity indices as columns. 
}
\details{
Let \emph{S_i} be the number of species in community \emph{i}, \eqn{n_{ij}}{n_ij} be the absolute abundance of species \emph{j} in community \emph{i}, \eqn{N_i} the sum of all species abundance in community \emph{i} (\eqn{N_i=\sum_j n_{ij}}{N_i=sum_j n_ij}; the sum of row \emph{i} in \code{comm}), \eqn{p_{ij}}{p_ij} the relative abundance of species \emph{j} in community \emph{i} (\eqn{p_{ij}=n_{ij}/N_i}{p_ij=n_ij/N_i}).
If \code{method="richness"}, the diversity index is the number of species.
If \code{method="GiniSimpson"}, the diversity index is that of Gini (1912) and Simpson (1949): \eqn{1-\sum_j p_{ij}^2}{1-sum_j p_ij^2}.
If \code{method="Simpson"}, the diversity index is (Simpson 1949): \eqn{1/\sum_j p_{ij}^2}{1/(sum_j p_ij^2)}.
If \code{method="Shannon"}, the diversity index is that of Shannon (1948) with neperian logarithm: \eqn{-\sum_j p_{ij}log(p_{ij})}{-sum_j p_ij log(p_ij)}.
If \code{method="Margalef"}, the diversity index is that of Margalef (1972): \eqn{(S_i-1)/log(N_i)}.
If \code{method="Menhinick"}, the diversity index is that of Menhinick (1964): \eqn{S_i/\sqrt{N_i}}{S_i/sqrt(N_i)}.
If \code{method="McIntosh"}, the diversity index is that of McIntosh (1967): \eqn{(N_i-\sqrt{\sum_j n_{ij}^2})/(N_i-\sqrt{N_i})}{(N_i-sqrt(sum_j n_ij^2))/(N_i-sqrt(N_i))}.
If one of the strings is "full", then all indices are calculated.
}

\references{
Gini, C. (1912) Variabilita e mutabilita. Studi economicoaguridici delle facoltta di giurizprudenza dell, Universite di Cagliari III, Parte II.

Magurran, A.E. (2004) Measuring biological diversity. Blackwell Publishing, Oxford, U.K.

Margalef, R. (1972) Homage to Evelyn Hutchinson, or why is there an upper limit to diversity? \emph{Transactions of the Connecticut Academy of Arts and Sciences}, \bold{44}, 211--235.

McIntosh, R.P. (1967) An index of diversity and the relation of certain concepts to diversity. \emph{Ecology}, \bold{48}, 392--404. 

Menhinick, E.F. (1964) A Comparison of some species-individuals diversity indices applied to samples of field insects. \emph{Ecology}, \bold{45}, 859--861.

Shannon, C.E. (1948) A mathematical theory of communication. \emph{Bell System technical journal}, \bold{27}, 379--423, 623--656.

Simpson, E.H. (1949) Measurement of diversity. \emph{Nature}, \bold{163}, 688.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
data(batcomm)
ab <- batcomm$ab
speciesdiv(ab)
}
\keyword{models}
