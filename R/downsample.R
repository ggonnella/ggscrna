#'
#' Downsample a Seurat object to a given number of cells
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

#' Downsample a Seurat object to a given number of cells
#'
#' @param so Seurat object
#' @param ncells Number of cells to downsample to
#' @param verbose Print verbose messages
#' @param tic Print timing information
#' @return Seurat object with the desired number of cells
#'         (or less if the original number of cells is less than ncells)
#' @export
#'
downsample_so <- function(so, ncells, verbose = TRUE, tic = TRUE) {
  if (tic) { tic("Downsampling Seurat object") }
  prev_n_cells <- length(so$cells)
  if (verbose) {
    print(paste0("Number of cells before downsampling: ",
                 prev_n_cells))
  }
  if (ncells >= prev_n_cells) {
    if (verbose) {
      print(paste0("No downsampling needed, as ",
                   ncells, " >= ", prev_n_cells))
    }
  }
  so <- subset(so, downsample = ncells)
  if (verbose) {
    print(paste0("Number of cells after downsampling: ",
                 length(so$cells)))
  }
  if (tic) { toc() }
  so
}
