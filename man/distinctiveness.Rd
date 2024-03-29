% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distinctiveness.R
\name{distinctiveness}
\alias{distinctiveness}
\title{Functional Distinctiveness on site-species matrix}
\usage{
distinctiveness(pres_matrix, dist_matrix, relative = FALSE)
}
\arguments{
\item{pres_matrix}{a site-species matrix (presence-absence or relative
abundances), with sites in rows and species in columns}

\item{dist_matrix}{a species functional distance matrix}

\item{relative}{a logical indicating if distinctiveness should be scaled
relatively to the community (scaled by max functional
distance among the species of the targeted community)}
}
\value{
a similar matrix from provided \code{pres_matrix} with Distinctiveness
values in lieu of presences or relative abundances, species absent from
communities will have an \code{NA} value (see \code{Note} section)
}
\description{
Computes functional distinctiveness from a site-species matrix (containing
presence-absence or relative abundances) of species with provided functional
distance matrix. The sites-species matrix should have \strong{sites} in
\strong{rows} and \strong{species} in \strong{columns}, similar to
\pkg{vegan} package defaults.
}
\details{
The Functional Distinctiveness of a species is the average functional
distance from a species to all the other in the given community. It is
computed as such:
\deqn{
   D_i = \frac{\sum_{j = 0, i \neq j}^N d_{ij}}{N-1},
   }{%
   D_i = ( \Sigma_(j = 0, i != j)^N d_ij) / (N - 1),
   }
with \eqn{D_i} the functional distinctiveness of species \eqn{i}, \eqn{N}
the total number of species in the community and \eqn{d_{ij}}{d_ij} the
functional distance between species \eqn{i} and species \eqn{j}.
\strong{IMPORTANT NOTE}: in order to get functional rarity indices between 0
and 1, the distance metric has to be scaled between 0 and 1.
}
\section{Note}{

Absent species should be coded by \code{0} or \code{NA} in input matrices.

When a species is alone in its community the functional distinctiveness
cannot be computed (denominator = 0 in formula), and its value is assigned
as \code{NaN}.

For speed and memory efficiency sparse matrices can be used as input of
the function using \code{as(pres_matrix, "dgCMatrix")} from the
\code{Matrix} package.
(see \code{vignette("sparse_matrices", package = "funrar")})
}

\examples{
data("aravo", package = "ade4")
# Site-species matrix
mat = as.matrix(aravo$spe)

# Compute relative abundances
mat = make_relative(mat)

# Example of trait table
tra = aravo$traits[, c("Height", "SLA", "N_mass")]
# Distance matrix
dist_mat = compute_dist_matrix(tra)

di = distinctiveness(pres_matrix = mat, dist_matrix = dist_mat)
di[1:5, 1:5]

# Compute distinctiveness for all species in the regional pool
# i.e., with all the species in all the communities
# Here considering each species present evenly in the regional pool
reg_pool = matrix(1, ncol = ncol(mat))
colnames(reg_pool) = colnames(mat)
row.names(reg_pool) = c("Regional_pool")

reg_di = distinctiveness(reg_pool, dist_mat)

}
