\name{specieseve}
\alias{specieseve}
\title{
Indices of Species Evenvess
}
\description{
The function \code{specieseve} calculates evenness indices that rely on relative or absolute species abundance. 
}
\usage{
specieseve(comm, method = "full", tol = 1e-8)
}
\arguments{
  \item{comm}{
a data frame or a matrix typically with communities as rows, species as columns and abundance as entry. 
}
  \item{method}{
a string or a vector of strings: one or several of "GiniSimpson", "Simpson", "Shannon", "Heip", "McIntosh", "SmithWilson", "full". See details.
}
  \item{tol}{
a tolerance threshold (a value between -\code{tol} and \code{tol} is considered equal to zero)}
}
\details{
Let \emph{S_i} be the number of species in community \emph{i}, \eqn{n_{ij}}{n_ij} be the absolute abundance of species \emph{j} in community \emph{i}, \eqn{N_i} the sum of all species abundance in community \emph{i} (\eqn{N_i=\sum_j n_{ij}}{N_i=sum_j n_ij}; the sum of row \emph{i} in \code{comm}), \eqn{p_{ij}}{p_ij} the relative abundance of species \emph{j} in community \emph{i} (\eqn{p_{ij}=n_{ij}/N_i}{p_ij=n_ij/N_i}).
If \code{method="GiniSimpson"}, the evenness index is that associated with Gini (1912) and Simpson (1949) diversity index: \eqn{(1-\sum_j p_{ij}^2)*S_i/(S_i-1)}{(1-\sum_j p_ij^2)*S_i/(S_i-1)}.
If \code{method="Simpson"}, the evenness index is (Simpson 1949; Magurran 2004): \eqn{(1/\sum_j p_{ij}^2)/S_i}{(1/sum_j p_ij^2)/S_i}.
If \code{method="Shannon"}, the evenness index is that associated with Shannon (1948) diversity index with neperian logarithm: \eqn{(-\sum_j p_{ij}ln(p_{ij}))/ln(S_i)}{(-sum_j p_ij ln(p_ij))/ln(S_i)}.
If \code{method="Heip"}, the evenness index is that of Heip (1974) (Magurran 2004): \eqn{[exp(-\sum_j p_{ij}log(p_{ij})) - 1]/(S_i-1)}{[exp(-sum_j p_ij log(p_ij)) - 1]/(S_i-1)}.
If \code{method="McIntosh"}, the evenness index is that of Pielou (1975) associated with McIntosh (1967) index of diversity: \eqn{(N_i-\sqrt{\sum_j n_{ij}^2})/(N_i-N_i/\sqrt{S_i})}{(N_i-sqrt(sum_j n_ij^2)/(N_i-N_i/sqrt(S_i))}.
If \code{method="SmithWilson"}, the Smith and Wilson (1996) evenness index is calculated (Magurran 2004): \eqn{1-[(2/\pi)(arctan(\sum_{j=1}^{S_i} (log n_{ij} - \sum_{k=1}^{S_i} log n_{ik}/S_i)^2/S_i))]}{1-[2/pi*arctan(sum_j (log(n_ij) - sum_k log(n_ik)/S_i)^2/S_i)]}. The function uses neperian logarithm for all indices. If one of the strings is "full", then all indices are calculated.
}
\value{
Function \code{specieseve} returns a matrix with communities as rows and the evenness indices as columns. 
}
\references{

Gini, C. (1912) \emph{Variabilita e mutabilita}. Studi economicoaguridici delle facoltta di giurizprudenza dell, Universite di Cagliari III, Parte II.

Heip, C. (1974) A new index measuring evenness. \emph{Journal of the Marine Biological Association UK}, \bold{54}, 555--557.

Magurran, A.E. (2004) \emph{Measuring biological diversity}. Oxford, UK: Blackwell Publishing.

McIntosh, R.P. (1967) An index of diversity and the relation of certain conepts to diversity. \emph{Ecology}, \bold{48}, 392--404. 

Pielou, E.C. (1975) \emph{Ecological diversity}. New York: Wiley InterScience.

Shannon, C.E. (1948) A mathematical theory of communication. \emph{Bell System technical journal}, \bold{27}, 379--423, 623--656.

Simpson, E.H. (1949) Measurement of diversity. \emph{Nature}, \bold{163}, 688.

Smith, B. and Wilson, J.B. (1996) A consumer's guide to evenness measures. \emph{Oikos}, \bold{76}, 70--82.

}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}
}
\examples{
data(batcomm)
ab <- batcomm$ab
specieseve(ab)
}
\keyword{models}
