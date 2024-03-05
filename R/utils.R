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

