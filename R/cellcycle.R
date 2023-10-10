#'
#' Cell cycle phase assignment
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

#'
#' Acknowledgements:
#' some of the code was originally derived from code written by Sebastien Mella
#'

#' Identify cell cycle phases
#'
#' @param so Seurat object
#' @param genes List of genes for cell cycle scoring
#' @param verbose Print verbose messages
#' @param tic Print timing information
#' @param plot Whether to plot the cell cycle phases using plot_cc_phase
#'
#' @return Seurat object with additional columns 'Phase', 'S.Score'
#'         and 'G2M.Score'
#'
id_cc_phase <- function(so, genes = Seurat::cc.genes.updated.2019,
                        verbose = TRUE, tic = TRUE, plot = FALSE) {
  if (tic) { tic("Identifying cell cycle phases") }
  s_genes = genes$s.genes
  g2m_genes = genes$g2m.genes
  s_genes_table <- as.data.frame(table(s_genes %in% rownames(so)))
  g2m_genes_table <- as.data.frame(table(g2m_genes %in% rownames(so)))
  so_norm <- Seurat::NormalizeData(so, assay = "RNA")
  so_norm <- Seurat::CellCycleScoring(so_norm,
                    s.features = s_genes[s_genes %in% rownames(so_norm)],
                    g2m.features = g2m_genes[g2m_genes %in% rownames(so_norm)])
  so$Phase <- as.factor(so_norm$Phase)
  so$S.Score <- so_norm$S.Score
  so$G2M.Score <- so_norm$G2M.Score
  if (verbose) { print(cc_phase_table(so)) }
  if (plot) { print(cc_phase_plot(so)) }
  if (tic) { toc() }
  so
}

#' Preamble to the plot/table functions
#'
#' @noRd
cc_phase_preamble <- function(so) {
  if (!"Phase" %in% colnames(so@meta.data)) {
    stop("Phase not found in metadata. Run id_cc_phase first.")
  }
  if (!is.null(sample)) {
    if (!sample %in% so$sample) {
      stop(paste0("Sample ", sample, " not found in Seurat object"))
    }
    so <- subset(so, sample == sample)
  }
  so
}

#' Compile a table of cell cycle phases
#'
#' @param so Seurat object
#' @param sample Sample name (optional) if the Seurat object contains
#'               multiple samples, only the specified sample is considered;
#'               default: consider all samples
#' @return Character vector with the table
#'
cc_phase_table <- function(so, sample = NULL) {
  so <- cc_phase_preamble(so)
  result <- character()
  result <- c(result, "Number of cells in each cell cycle phase:")
  result <- c(result, "Phase\tNumber\tPercentage")
  for (phase in levels(so$Phase)) {
    perc = sum(so$Phase == phase) / length(so$Phase) * 100
    result <- c(result, paste0(phase, "\t", sum(so$Phase == phase), "\t",
                               round(perc, 1), "%"))
  }
  result
}

#' Create a plot of cell cycle phases
#'
#' @param so Seurat object
#' @param sample Sample name (optional) if the Seurat object contains
#'               multiple samples, only the specified sample is considered;
#'               default: consider all samples
#' @return ggplot2 object
#'
cc_phase_plot <- function(so, sample = NULL) {
  so <- cc_phase_preamble(so)
  ccp_df <- as.data.frame.array(table(so$Phase))
  ccp_df$phase <- rownames(ccp_df)
  so <- Seurat::FindVariableFeatures(so, selection.method = "vst",
                                     array = "RNA", nfeatures = 2000,
                                     verbose = FALSE)
  so <- Seurat::ScaleData(so, assay = "RNA", verbose = FALSE)
  so <- Seurat::RunPCA(so, assay = "RNA", reduction.name = "pcaRNA",
               reduction.key = "pcaRNA_", verbose = FALSE)
  ccg_colors <-
    RColorBrewer::brewer.pal(length(levels(so$Phase)), name = "Set1")
  plt1 <- Seurat::DimPlot(so, reduction = "pcaRNA", group.by= "Phase",
                         cols = ccg_colors)
  plt2 <- Seurat::DimPlot(so, reduction = "pcaRNA", group.by= "Phase",
                         split.by = "Phase", cols = ccg_colors)
  plt <- gridExtra::grid.arrange(plt1, plt2, ncol = 2)
  plt
}

