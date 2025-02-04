% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{consolidate_matrix_rows}
\alias{consolidate_matrix_rows}
\title{Consolidate rows of a matrix}
\usage{
consolidate_matrix_rows(matrix, rows_to_sum, consolidated_label)
}
\arguments{
\item{matrix}{Matrix. The matrix to be consolidated}

\item{rows_to_sum}{Vector of strings. The row names to be summed}

\item{consolidated_label}{String. The label of the consolidated row}
}
\value{
Matrix with the specified rows summed and the
                           consolidated row added at the end. The rows in
                           rows_to_sum are removed from the matrix.
}
\description{
Consolidate rows of a matrix
}
