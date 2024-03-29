% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stack_matrix.R
\name{stack_to_matrix}
\alias{stack_to_matrix}
\alias{tidy_to_matrix}
\title{Stacked (= tidy) data.frame to matrix}
\usage{
stack_to_matrix(
  my_df,
  col_to_row,
  col_to_col,
  col_value = NULL,
  sparse = FALSE
)
}
\arguments{
\item{my_df}{data.frame you want to transform in matrix}

\item{col_to_row}{character vector of the name of the data.frame column you
want to put into matrix rows}

\item{col_to_col}{character vector of the name of the data.frame column you
want to be as columns in matrix}

\item{col_value}{(optional, default = \code{NULL}) character vector indicating
the name of a column coding the values that will be put in
the matrix}

\item{sparse}{(optional, default = \code{FALSE}) logical indicating whether to
return a sparse matrix (if \code{TRUE} requires
\href{https://cran.r-project.org/package=tidytext}{\code{tidytext}}
package)}
}
\value{
a matrix with given \code{col_to_row} column in rows and \code{col_to_col}
column in columns. If some cells are not present in the data.frame (e.g. some
species not present at some sites), the matrix will have a \code{NA} value.
}
\description{
Passes from a stacked (= tidy) data.frame to a matrix. \code{tidy_to_matrix()} is
an alias of this function.
}
\examples{
example = data.frame("sites" = c(rep("1", 3), rep("2", 2)),
 "species" = c("A", "B", "C", "B", "D"),
  "abundance" = c(0.33, 0.33, 0.33, 0.4, 0.6))

mat = stack_to_matrix(example, "sites", "species", "abundance")
mat

}
\seealso{
\code{\link[=matrix_to_stack]{matrix_to_stack()}} for the reverse operation
}
