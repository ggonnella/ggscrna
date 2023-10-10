#'
#' Identify gene sets
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

geneset_rx <- list()
geneset_rx$TCR  <- "^TR[AB][VC]"
geneset_rx$BCR  <- "^IG[HKL][VC]"
geneset_rx$MT   <- "^MT-"
geneset_rx$RIBO <- "^M?RP.*"

#' Identify genes in a gene set
#'
#' @param so a Seurat object
#' @param geneset_name a string, one of "TCR", "BCR", "MT", "RIBO"
#' @param verbose a boolean, if TRUE print the number of genes in the set
#' @param tic a boolean, if TRUE print the time taken to identify the genes
#' @return a vector of gene names
#'
id_geneset <- function(so, geneset_name, verbose=TRUE, tic=FALSE) {
  if (tic) { tic(paste("Identifying", geneset_name, "genes")) }
  if (!is(so, "Seurat")) {
    stop("so must be a Seurat object")
  }
  if (!geneset_name %in% names(geneset_rx)) {
    stop(paste("geneset_name must be one of", names(geneset_rx)))
  }
  rx <- geneset_rx[[geneset_name]]
  geneset <- grep(rx, rownames(so), value = TRUE)
  if (verbose) {
    print(paste("Number of", geneset_name, "genes:", length(geneset)))
  }
  if (tic) { toc() }
  geneset
}

#' Remove gene set from a Seurat object
#'
#' @param so a Seurat object
#' @param geneset_name a string, one of "TCR", "BCR", "MT", "RIBO"
#' @param verbose a boolean, if TRUE print the number of genes in the set
#' @param tic a boolean, if TRUE print the time taken to remove the genes
#' @return a Seurat object
#'
rm_geneset <- function(so, geneset_name, verbose=TRUE, tic=FALSE) {
  if (tic) { tic(paste("Identifying", geneset_name, "genes")) }
  if (!is(so, "Seurat")) {
    stop("so must be a Seurat object")
  }
  geneset <- id_geneset(so, geneset_name, verbose=verbose, tic=FALSE)
  if (length(geneset) > 0) {
    so <- so[!rownames(so) %in% geneset,]
    if (verbose) {
      print(paste("Removed", length(geneset), geneset_name, "genes"))
    }
  } else {
    if (verbose) {
      print(paste("No", geneset_name, "genes found"))
    }
  }
  if (tic) { toc() }
  so
}

#' Remove multiple genesets based on a list of parameters
#'
#' @param so a Seurat object
#' @param params parameters, either a list or environment.
#            if an element is one of "TCR", "BCR", "MT",
#'           "RIBO" and the value is 'rm', then the corresponding gene set is
#'           removed
rm_genesets <- function(so, params, verbose=TRUE, tic=TRUE) {
  if (tic) { tic("Removing genesets") }
  if (!is(so, "Seurat")) {
    stop("so must be a Seurat object")
  }
  for (geneset_name in names(geneset_rx)) {
    if (!geneset_name %in% names(params)) {
      next
    }
    if (params[[geneset_name]] == "rm") {
      so <- rm_geneset(so, geneset_name, verbose=verbose, tic=FALSE)
    }
  }
  if (tic) { toc() }
  so
}
