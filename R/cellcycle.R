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
#'
#' @return Seurat object with additional columns 'Phase', 'S.Score'
#'         and 'G2M.Score'
#'
id_cc_phase <- function(so, genes = Seurat::cc.genes.updated.2019) {
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
  so
}

#' Plot cell cycle phases
#'
plot_cc_phase <- function(so, sample = NULL) {
  if (!is.null(sample)) {
    so <- subset(so, sample == sample)
  }
  so$Phase = as.factor(so$Phase)
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
  plt1 = Seurat::DimPlot(so, reduction = "pcaRNA", group.by= "Phase",
                         cols = ccg_colors)
  plt2 = Seurat::DimPlot(so, reduction = "pcaRNA", group.by= "Phase",
                         split.by = "Phase", cols = ccg_colors)
  gridExtra::grid.arrange(plt1, plt2, ncol = 2)
}

