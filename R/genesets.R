#'
#' Identify gene sets
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023-2024

genesets <- c("TCR", "BCR", "MT", "RIBO")

regex_geneset <- list()
regex_geneset$TCR  <- "^TR[AB][VC]"
regex_geneset$BCR  <- "^IG[HKL][VC]"
regex_geneset$MT   <- "^MT-"
regex_geneset$RIBO <- "^M?RP.*"

#' Identify genes in a gene set
#'
#' @param so            the Seurat object
#' @param geneset_name  String. Geneset name
#'                      it must be one of the elements of 'genesets'
#' @param verbose       Boolean. TRUE: be verbose
#' @param tic           Boolean. TRUE: show the running time
#' 
#' @return              Vector. Gene names
#'
id_geneset <- function(so, geneset_name, verbose=TRUE, tic=FALSE) {
  if (tic)
    tic(paste("Identifying", geneset_name, "genes"))

  if (!is(so, "Seurat"))
    stop("so must be a Seurat object")

  if (!geneset_name %in% genesets)
    stop(paste("geneset_name must be one of", genesets))

  genes <- grep(regex_geneset[[geneset_name]], rownames(so), value = TRUE)

  if (verbose)
    print(paste("Number of", genes_name, "genes:", length(genes)))

  if (tic)
    toc()

  genes
}

#' Remove gene set from a Seurat object
#'
#' @param so            the Seurat object
#' @param geneset_name  String. Geneset name
#'                      it must be one of the elements of 'genesets'
#' @param verbose       Boolean. TRUE: be verbose
#' @param tic           Boolean. TRUE: show the running time
#' 
#' @return              the Seurat object without the genes in the specified geneset
#'
rm_geneset <- function(so, geneset_name, verbose=TRUE, tic=FALSE) {
  if (tic)
    tic(paste("Remove", geneset_name, "genes"))

  if (!is(so, "Seurat"))
    stop("so must be a Seurat object")

  if (!geneset_name %in% genesets)
    stop(paste("geneset_name must be one of", genesets))

  genes <- id_genes(so, geneset_name, verbose=verbose, tic=FALSE)
  if (length(genes) > 0) {
    so <- so[!rownames(so) %in% genes,]
    if (verbose)
      print(paste("Removed", length(genes), geneset_name, "genes"))
  } else {
    if (verbose)
      print(paste("No", geneset_name, "genes found"))
  }

  if (tic)
    toc()

  so
}

#' Remove multiple genesets based on a list of parameters
#' 
#' It calls rm_geneset for each geneset in the parameters list or env,
#' where a parameter is called after a geneset and its value is 'rm'
#'
#' @param so       the Seurat object
#' @param params   List or Environment. Parameters
#' @param verbose  Boolean. TRUE: be verbose
#' @param tic      Boolean. TRUE: show the running time
#' 
#' @return         the Seurat object without the genes in the specified genesets
#'
rm_genesets <- function(so, params, verbose=TRUE, tic=TRUE) {
  if (tic)
    tic("Removing genesets")

  if (!is(so, "Seurat"))
    stop("so must be a Seurat object")

  for (geneset_name in genesets) {
    if (!geneset_name %in% names(params))
      next
    if (params[[geneset_name]] == "rm")
      so <- rm_geneset(so, geneset_name, verbose=verbose, tic=FALSE)
  }

  if (tic)
    toc()

  so
}
