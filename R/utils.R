#'
#' Helper functions for the giosc package
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023


#' Include a Rmd file in another
#'
#' @param file Filename of the file to be included
#'
source_rmd <- function(file, ...) {
  tmp_file <- tempfile(fileext = ".R")
  on.exit(unlink(tmp_file), add = TRUE)
  knitr::purl(file, output = tmp_file)
  source(file <- tmp_file, ...)
}

