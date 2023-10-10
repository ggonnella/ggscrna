#'
#' Cell type assignment
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

#'
#' Acknowledgements:
#' some of the code was originally derived from code written by Sebastien Mella
#'

#' Cell type assignment using SingleR
#'
#' @param so Seurat object
#'
#' @return Seurat object with additional columns 'SingleR.best_hit',
#'         'SingleR.best_score' and 'SingleR.best_match'
#'
assign_cell_types <- function(so) {
  ref <- celldex::MonacoImmuneData()
  sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = so@assays$RNA@counts))
  sce <- scater::logNormCounts(sce)
  for (predtype in c("fine", "main")) {
    if (predtype == "fine") {
      lbls <- ref$label.fine
    } else {
      lbls <- ref$label.main
    }
    pred <- SingleR::SingleR(test = sce, ref = ref, labels = lbls,
                    assay.type.test = "logcounts")
    pred_tab <- as.data.frame(pred$scores[, colnames(pred$scores) %in%
                  unique(pred$labels)])
    pred_tab$pruned.labels <- pred$pruned.labels
    colnames(pred_tab) <- paste(colnames(pred_tab),
                                   "singleR", predtype, sep = "_")
    rownames(pred_tab) <- rownames(pred)
    non_pred_cols <- setdiff(colnames(so@meta.data), colnames(pred_tab))
    so@meta.data <- so@meta.data[, non_pred_cols]
    so@meta.data <- cbind(so@meta.data, pred_tab)
    lbl <- paste0("pruned.labels_singleR_", predtype)
    so[[lbl]][is.na(so[[lbl]])] <- "unchar"
    lev <- names(sort(table(so[[lbl]]), decreasing = TRUE))
    lev <- c(lev[-which(lev == "unchar")], "unchar")
    so_vector <- as.vector(so[[lbl]][,1])
    so@meta.data[[lbl]] <- factor(so_vector, levels = lev)
  }
  so
}
