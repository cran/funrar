% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funrar.R
\name{funrar}
\alias{funrar}
\title{Compute all Functional Rarity Indices from Matrices}
\usage{
funrar(pres_matrix, dist_matrix, rel_abund = FALSE)
}
\arguments{
\item{pres_matrix}{a site-species matrix (presence-absence or relative
abundances), with sites in rows and species in columns}

\item{dist_matrix}{a species functional distance matrix}

\item{rel_abund}{logical (\code{TRUE} or \code{FALSE}) indicating if
site-species matrix contain relative abundances values or only
presence-absence data (default = \code{FALSE})}
}
\value{
A list of 3 objects (or 4 if \code{rel_abund = TRUE}):
\describe{
\item{\strong{Ui}}{a vector containing uniqueness values per species,}
\item{\strong{Di}}{a site-species matrix with functional distinctiveness
values per species per site,}
\item{\strong{Ri}}{a vector containing geographical restrictedness values
per species,}
}
and if \code{rel_abund = TRUE},
\describe{
\item{\strong{Si}}{a site-species matrix with scarcity values per species
per site.}
}
}
\description{
From a site-species matrix and functional distance matrix compute all indices
included in the package: functional uniqueness (regional, functional),
functional distinctiveness (local, functional), geographical restrictedness
(regional, extent), scarcity (local, abundance). \strong{Note}: scarcity can
only be computed if relative abundances are provided in the site-species
matrix.
}
\seealso{
\code{\link[=uniqueness]{uniqueness()}}, \code{\link[=distinctiveness]{distinctiveness()}}, \code{\link[=restrictedness]{restrictedness()}},
\code{\link[=scarcity]{scarcity()}}
}
