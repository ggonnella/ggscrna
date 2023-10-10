#'
#' Cell type assignment
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

#'
#' Acknowledgements:
#' some of the code was originally derived from code written by Sebastien Mella
#'

#' Get label set from a celldex ref dataset using a string identifier
#' @noRd
celldex_get_labels <- function(ref, lbl_id) {
  if (lbl_id == "main") {
    lbls <- ref$label.main
  } else if (lbl_id == "fine") {
    lbls <- ref$label.fine
  } else if (lbl_id == "ont") {
    lbls <- ref$label.ont
  } else {
    stop("lbl_id must be one of: main, fine, ont")
  }
  lbls
}

prepare_colnames <- function(cns, predtype, refname) {
  cns <- paste(cns, "singleR", predtype, refname, sep = "_")
  cns <- gsub(" ", "_", cns)
  cns
}

#' Transfer SingleR results to Seurat object
#'
#' pred contains:
#'   scores, tuning.scores
#'   first.labels, labels, pruned.labels
#'
#' this function:
#' (1) removes any previous results
#' (3) transfers the pruned.labels
#' (2) transfers scores used in pruned.labels
#'
#' columns are prefixed by singleR, predtype, refname
#' spaces in the column names are replaced by underscores
#'
#' @noRd
singleR_res_to_so_metadata <- function(so, pred) {
  assigned_labels <- unique(pred$pruned.labels)
  is_score_col_used <- colnames(pred$scores) %in% assigned_labels
  pred_tab <- as.data.frame(pred$scores[,is_score_col_used])
  pred_tab$pruned.labels <- pred$pruned.labels
  rownames(pred_tab) <- rownames(pred)
  colnames(pred_tab) <- prepare_colnames(colnames(pred_tab), predtype, refname)
  # remove previous results if any
  unused_cols <- colnames(pred$scores)[!is_score_col_used]
  unused_cols <- prepare_colnames(unused_cols, predtype, refname)
  all_cols <- c(colnames(pred_tab), unused_cols)
  non_pred_cols <- setdiff(colnames(so@meta.data), all_cols)
  so@meta.data <- so@meta.data[, non_pred_cols]
  # add new results
  so@meta.data <- cbind(so@meta.data, pred_tab)
  so
}

#' Postprocess SingleR labels
#'   convert labels into factors;
#'   convert NA labels to 'Uncharacterized'
#' @noRd
postprocess_singleR_labels <- function(so, col, ulbl = "Uncharacterized") {
  so[[col]][is.na(so[[col]])] <- ulbl
  lev <- names(sort(table(so[[col]]), decreasing = TRUE))
  lev <- c(lev[-which(lev == ulbl)], ulbl)
  so_vector <- as.vector(so[[col]][,1])
  so@meta.data[[col]] <- factor(so_vector, levels = lev)
  so
}


#' Cell type assignment using SingleR and celldex reference
#'
#' @param so Seurat object
#' @param refs. Reference datasets of celldex to use.
#'          Default: MonacoImmuneData, ImmGenData
#'          Other generic datasets:
#'            HumanPrimaryCellAtlasData
#'            MouseRNAseqData
#'          Other immunology datasets:
#'            BlueprintEncodeData:
#'               low resolution
#'               best option for: fast interpretation of results
#'            ImmGenData:
#'               very high resolution
#'            DatabaseImmuneCellExpressionData:
#'               good resolution: CD4+ T cells
#'               low resolution: B cells
#'               lacking: CD4+ central memory, effector memory, dendritic cells
#'            NovershternHematopoieticData:
#'              good resolution: myeloid cells, progenitors, NK, erythroid,
#'                               granulocytes
#'              best option for: bone marrow samples
#' @param refnames. Names of the reference datasets.
#'                 Default: "monaco", "immgen"
#' @param predtypes. Types of labels to use.
#'                   Default: "fine", "main", "ont"
#' @param getlbl. Function to get the labels from the reference object.
#'                using fun(ref, predtype). Defaults to an internal function
#'                which can be used with celldex datasets.
#' @return Seurat object with additional meta.data columns:
#'           pruned.labels_singleR_<refname>_<predtype>
#'
assign_cell_types <- function(so,
                              refs = c(celldex::MonacoImmuneData()),
                              refnames = c("monaco"),
                              predtypes = c("fine", "main", "ont"),
                              getlbl = celldex_get_labels) {
  if (length(ref) != length(refnames)) {
    stop("length(ref) != length(refnames)")
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = so@assays$RNA@counts))
  sce <- scater::logNormCounts(sce)
  for (i in 1:length(ref)) {
    ref <- refs[[i]]
    refname <- refnames[[i]]
    for (predtype in predtypes) {
      pred <- SingleR::SingleR(test = sce, ref = ref,
                               labels = getlbl(ref, predtype),
                               assay.type.test = "logcounts")
      so <- singleR_res_to_so_metadata(so, pred)
      col <- paste0("pruned.labels_singleR_", predtype, "_", refname)
      postprocess_singleR_labels(so, col)
    }
  }
  so
}
