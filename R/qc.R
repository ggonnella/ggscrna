#'
#' Quality control functions.
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

#'
#' Acknowledgements:
#' some of the code was originally derived from code written by Sebastien Mella
#'

# define %>%
`%>%` <- dplyr::`%>%`

#' Compute additional metrics on the Seurat object
#'
#' The metrics are those explained in:
#' https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
#'
#' ## Notes
#' - mt genes computation by MT- prefix assumes this is a human dataset
#'
#' @param so SeuratObject for which to compute the metrics
#'
#' @return The SeuratObject, with additional metrics in ``@meta.data``
#"         and renaming some of the default metrics.
prepare_QC <- function(so) {
  so$mitoRatio <- Seurat::PercentageFeatureSet(object = so, pattern = "^MT-")
  so$mitoRatio <- so@meta.data$mitoRatio / 100
  so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
  so$cells <- rownames(so@meta.data)
  so
}

#' Add QC information from scater library to Seurat object
#'
#' @param so Seurat object
#'
#' @return Seurat object with additional QC metrics
scater_qc <- function(so) {
  sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = so@assays$RNA@counts))
  is_mito <- grep("^MT-", rownames(sce))
  qcstats <- scater::perCellQCMetrics(sce, subsets=list(Mito=is_mito))
  filt_ad <- data.frame(
    qc.lib.low = scater::isOutlier(qcstats$sum, log = TRUE, type = "lower"),
    qc.lib.high = scater::isOutlier(qcstats$sum,
                            log = FALSE, type = "higher", nmads = 4),
    qc.nexpr.low = scater::isOutlier(qcstats$detected, log = TRUE,
                                     type = "lower"),
    qc.nexpr.high = scater::isOutlier(qcstats$detected,
                              log = FALSE, type = "higher", nmads = 4),
    qc.mito.high = scater::isOutlier(qcstats$subsets_Mito_percent,
                             log = FALSE, type = "higher")
  )
  qcstats <- cbind(qcstats, filt_ad)
  names(qcstats) <- paste0("scater_", names(qcstats))
  so@meta.data <- cbind(so@meta.data, qcstats)
  withr::local_seed(1234)
  stats <- cbind(log10(so$nCount_RNA), log10(so$nFeature_RNA), so$mitoRatio)
  outlying <- robustbase::adjOutlyingness(stats, only.outlyingness = TRUE)
  multi.outlier <- scater::isOutlier(outlying, type = "higher")
  so@meta.data <- cbind(so@meta.data, multi.outlier)
  so
}

#' QC metric distribution histogram
#'
#' An histogram plot of the distribution
#' for one of the QC metrics in a Seurat object
#'
#' @param so Seurat object
#' @param measure String. The name of the metric to plot
#'
#' @return A ggplot object
#'
qc_distri_plot <- function(so, measure) {
  ggplot2::ggplot(so@meta.data, ggplot2::aes_string(x=measure)) +
    ggplot2::geom_histogram(bins=100) +
  ggplot2::ggtitle(paste(measure, "Distribution")) +
  ggplot2::xlab(measure) + ggplot2::ylab("Frequency")
}


#' Plot number of cells in different samples
#'
#' @param so Seurat object
#'
histogram_n_cells <- function(so) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(x=sample, fill=sample)) +
    	ggplot2::geom_bar() +
    	ggplot2::theme_classic() +
    	ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
    	ggplot2::theme(plot.title =
                     ggplot2::element_text(hjust=0.5, face="bold")) +
    	ggplot2::ggtitle("Number of cells")
}

#' Plot number of UMIs per cell in different samples
#'
#' @param so Seurat object
#'
density_plot_n_umis <- function(so, plot_nrows=0) {
  so@meta.data %>%
      ggplot2::ggplot(ggplot2::aes(color=sample, x=nCount_RNA, fill=sample)) +
    	ggplot2::geom_density(alpha = 0.2) +
    	ggplot2::scale_x_log10() +
    	ggplot2::theme_classic() +
    	ggplot2::ylab("cell density") +
    	ggplot2::geom_vline(xintercept = 100) +
    	ggplot2::geom_vline(xintercept = 500) +
    	ggplot2::geom_vline(xintercept = 1000) +
      ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
    	ggplot2::ggtitle("Number of UMIs per cell")
}

#' Boxplot of the number of UMIs per cell
#'
#' @param so Seurat object
#'
boxplot_n_umis <- function(so) {
  so@meta.data %>%
    	ggplot2::ggplot(
          ggplot2::aes(x=sample, y=log10(nCount_RNA), fill=sample)) +
    	ggplot2::geom_boxplot() +
    	ggplot2::theme_classic() +
    	ggplot2::theme(axis.text.x =
          ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
    	ggplot2::theme(plot.title =
          ggplot2::element_text(hjust=0.5, face="bold")) +
    	ggplot2::ggtitle("Number of UMIs per cell")
}

#' Plot number of genes per cell in different samples
#'
#' @param so Seurat object
#'
density_plot_n_genes <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(color=sample, x=nFeature_RNA, fill=sample)) +
    	ggplot2::geom_density(alpha = 0.2) +
    	ggplot2::theme_classic() +
    	ggplot2::scale_x_log10() +
    	ggplot2::geom_vline(xintercept = 300) +
      ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
    	ggplot2::ggtitle("Number of genes per cell")
}


#' Boxplot of the number of genes per cell
#'
#' @param so Seurat object
#'
boxplot_n_genes <- function(so) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(x=sample,
                                   y=log10(nFeature_RNA), fill=sample)) +
    	ggplot2::geom_boxplot() +
    	ggplot2::theme_classic() +
    	ggplot2::theme(axis.text.x =
        ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
    	ggplot2::theme(plot.title =
        ggplot2::element_text(hjust=0.5, face="bold")) +
    	ggplot2::ggtitle("Number of genes per cell")
}

dotplot_n_umis_genes_mito <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(x=nCount_RNA,
                                   y=nFeature_RNA, color=mitoRatio)) +
    	ggplot2::geom_point() +
  	  ggplot2::scale_colour_gradient(low = "gray90", high = "black") +
    	ggplot2::stat_smooth(method=lm) +
    	ggplot2::scale_x_log10() +
    	ggplot2::scale_y_log10() +
    	ggplot2::theme_classic() +
    	ggplot2::geom_vline(xintercept = 500) +
    	ggplot2::geom_hline(yintercept = 250) +
    	ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
      ggplot2::ggtitle("N. UMIs vs N. genes and mit.genes ratio")
}

density_plot_mito_ratio <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot2::ggplot(ggplot2::aes(color=sample, x=mitoRatio, fill=sample)) +
    	ggplot2::geom_density(alpha = 0.2) +
    	ggplot2::scale_x_log10() +
    	ggplot2::theme_classic() +
    	ggplot2::geom_vline(xintercept = 0.2) +
      ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
      ggplot2::ggtitle("Mitochondrial gene ratio per cell")
}

density_plot_complexity <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot2::ggplot(
        ggplot2::aes(x=log10GenesPerUMI, color=sample, fill=sample)) +
    	ggplot2::geom_density(alpha = 0.2) +
    	ggplot2::theme_classic() +
    	ggplot2::geom_vline(xintercept = 0.8) +
      ggplot2::facet_wrap(~sample, nrow=plot_nrows) +
      ggplot2::ggtitle("Transcriptional complexity")
}

#' Create and print several QC metrics plots
#'
#' Creates several QC metrics plots for a Seurat object
#' and prints them to the screen
#'
#' @param so Seurat object
#'
show_qc_plots <- function(so, plot_nrows=0) {
  print(histogram_n_cells(so))
  print(density_plot_n_umis(so, plot_nrows))
  print(boxplot_n_umis(so))
  print(density_plot_n_genes(so, plot_nrows))
  print(boxplot_n_genes(so))
  print(density_plot_mito_ratio(so, plot_nrows))
  print(dotplot_n_umis_genes_mito(so, plot_nrows))
  print(density_plot_complexity(so, plot_nrows))
}

#' Plot the highest expressed genes
#'
#' @param so Seurat object
#' @param n_genes Integer. Number of genes to plot
#' @param sample String. Keep only specified sample. Default: keep all.
#'
#' @return A ggplot object
h_expr_genes_plot <- function(so, n_genes = 10, sample = NULL) {
  title <- "Highest expressed genes"
  count_matrix <- as.matrix(so@assays$RNA@counts)
  if (!is.null(sample)) {
    count_matrix <- count_matrix[, grepl(sample, colnames(count_matrix))]
    title <- paste0(title, " (Sample ", sample, ")")
  }
  values <- head(order(matrixStats::rowSums2(count_matrix),
                       decreasing = TRUE), n_genes)
  sub_mat <- as.data.frame(count_matrix[values, ])
  sub_mat$feature <- factor(rownames(sub_mat), levels = rev(rownames(sub_mat)))
  sub_mat_long_fmt <- melt(sub_mat, id.vars = "feature")
  plot_colors <- grDevices::colorRampPalette(
                   rev(ggsci::pal_futurama("planetexpress")(12)))(n_genes)
  plt <- ggplot2::ggplot(sub_mat_long_fmt,
           ggplot2::aes(x = value, y = feature)) +
         ggplot2::geom_boxplot(aes(fill = feature), alpha = .75) +
         ggplot2::scale_fill_manual(values = plot_colors) +
         ggplot2::theme_light() +
         ggplot2::guides(fill="none") +
         ggplot2::xlab("Count") +
         ggplot2::ylab("Gene") +
         ggplot2::ggtitle(title)
  plt
}
