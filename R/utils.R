#'
#' Helper functions for the giosc package
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023-2024

#' Call the code of a Rmd file in another Rmd or a R script.
#' The code of the Rmd file is extracted as R script and then sourced. 
#'
#' @param file Filename of the file to be included
#'
source_rmd <- function(file, ...) {
  tmp_file <- tempfile(fileext = ".R")
  on.exit(unlink(tmp_file), add = TRUE)
  knitr::purl(file, output = tmp_file)
  source(file <- tmp_file, ...)
}

#' Get the path to a step directory
#'
#' @param steps_dir  String.  Path to the "steps" directory.
#' @param step_n     String.  Prefix of the step.
#'
#' @return           String.  Path to the step directory.
#' 
get_step_dir <- function(steps_dir, step_n) {
  all_files      <- list.files(path = steps_dir)
  matching_files <- all_files[grepl(paste0("^", step_n), all_files)]

  # there should be only one match
  if (length(matching_files) == 0)
    stop(paste0("No entry matches '",       step_n, "' in '", steps_dir, "'"))
  else if (length(matching_files) > 1)
    stop(paste0("Multiple entries match '", step_n, "' in '", steps_dir, "'"))

  paste0(steps_dir, "/", matching_files[1], "/")
}

#' Consolidate rows of a matrix
#' 
#' @param matrix              Matrix. The matrix to be consolidated
#' @param rows_to_sum         Vector of strings. The row names to be summed
#' @param consolidated_label  String. The label of the consolidated row
#' 
#' @return                    Matrix with the specified rows summed and the
#'                            consolidated row added at the end. The rows in
#'                            rows_to_sum are removed from the matrix.
#' 
consolidate_matrix_rows <- function(matrix, rows_to_sum, consolidated_label) {
  rows_to_sum <- rows_to_sum[rows_to_sum %in% rownames(matrix)]
  if (length(rows_to_sum) == 0)
    return(matrix)

  summed_rows <- colSums(as.matrix(matrix[rows_to_sum, ], drop=FALSE))
  new_matrix <- matrix[!rownames(matrix) %in% rows_to_sum, ]
  new_matrix <- rbind(new_matrix, summed_rows)
  rownames(new_matrix)[nrow(new_matrix)] <- consolidated_label

  new_matrix
}
