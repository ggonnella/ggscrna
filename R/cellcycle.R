#'
#' Cell cycle phase assignment
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023-2024

#' Identify cell cycle phases
#'
#' @param so       Seurat object
#' @param genes    List of genes for cell cycle scoring
#' @param verbose  Print verbose messages
#' @param tic      Print timing information
#' @param plot     Whether to plot the cell cycle phases using plot_cc_phase
#'
#' @return         Seurat object with additional metadata:
#'                 'Phase', 'S.Score', 'G2M.Score'
#'
id_cc_phase <- function(so, genes = Seurat::cc.genes.updated.2019,
                        verbose = TRUE, tic = TRUE, plot = FALSE) {
  if (tic) { tic("Identifying cell cycle phases") }

  so_copy <- Seurat::NormalizeData(so, assay = "RNA")
  so_copy <- Seurat::CellCycleScoring(so_copy,
                    s.features =   genes$s.genes[  genes$s.genes   %in% rownames(so_copy)],
                    g2m.features = genes$g2m.genes[genes$g2m.genes %in% rownames(so_copy)])

  # add metadata to the original object
  so$Phase     <- as.factor(so_copy$Phase)
  so$S.Score   <- so_copy$S.Score
  so$G2M.Score <- so_copy$G2M.Score

  if (verbose) { print(cc_phase_table(so)) }
  if (plot)    { print(cc_phase_plot(so))  }
  if (tic)     { toc()                     }

  so
}

#' Check that the phase metadata is available
#'
#' @noRd
chkavail_phase <- function(so) {
  if (!"Phase" %in% colnames(so@meta.data)) {
    stop("Phase not found in metadata. Run id_cc_phase first.")
  }
}

#' Check that a sample exists and extract it
#'
extract_subsample <- function(so, sample) {
    if (!sample %in% so$sample)
      stop(paste0("Sample ", sample, " not found in Seurat object"))
    subset(so, sample == sample)
}

#' Compile a table of cell cycle phases
#'
#' @param so      Seurat object
#' @param sample  Sample name (optional)
#'                if provided: process only specified sample
#'                default: consider all samples
#' 
#' @return        Character vector with the table
#'
cc_phase_table <- function(so, sample = NULL) {
  chkavail_phase(so)
  if (!is.null(sample))
    so <- extract_subsample(so, sample)

  result <- character()
  result <- c(result, "# Number of cells in each cell cycle phase")
  result <- c(result, "Phase\tNumber\tPercentage")

  for (phase in levels(so$Phase)) {
    total = sum(so$Phase == phase)
    perc = round(sum(so$Phase == phase) / length(so$Phase) * 100, 1)
    perc = paste0(perc, "%")
    result <- c(result, paste0(phase, "\t", total, "\t", perc))
  }

  result
}

#' Create a plot of cell cycle phases
#'
#' @param so       Seurat object
#' @param sample   Sample name (optional)
#'                 if provided: process only specified sample
#'                 default: consider all samples
#' @param outfile  Output file (optional)
#'                 if provided: save the plot to the file
#'                 default: return the plot, but do not save it to file
#' 
#' @return         ggplot2 object
#'
cc_phase_plot <- function(so, sample = NULL, outfile = NULL) {
  chkavail_phase(so)
  if (!is.null(sample))
    so <- extract_subsample(so, sample)

  # preprocess and create PCA
  so_copy <- Seurat::FindVariableFeatures(so, selection.method = "vst",
                                          array = "RNA", nfeatures = 2000,
                                          verbose = FALSE)
  so_copy <- Seurat::ScaleData(so_copy, assay = "RNA", verbose = FALSE)
  so_copy <- Seurat::RunPCA(so_copy, assay = "RNA",
                            reduction.name = "pcaRNA",
                            reduction.key = "pcaRNA_",
                            verbose = FALSE)

  # create palette
  n_phases <-   length(levels(so_copy$Phase))
  ccg_colors <- RColorBrewer::brewer.pal(n_phases, name = "Set1")

  # create a double plot, (1) all phases and (2) split by phases
  plt1 <- Seurat::DimPlot(so_copy, reduction = "pcaRNA", group.by= "Phase",
                         cols = ccg_colors)
  plt2 <- Seurat::DimPlot(so_copy, reduction = "pcaRNA", group.by= "Phase",
                         split.by = "Phase", cols = ccg_colors)
  plt <- gridExtra::grid.arrange(plt1, plt2, ncol = 2)

  # save to file if outfile is provided
  if (!is.null(outfile))
    ggsave(outfile, plt, width = 10, height = 5)

  plt
}

