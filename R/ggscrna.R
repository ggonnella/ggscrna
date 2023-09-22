#'
#' Helper functions for scRNA-seq analysis with Seurat
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

library(Seurat)
library(dplyr)
library(ggplot2)
library(glue)
library(knitr)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(matrixStats)
library(scater)

#' Get sample IDs from a sample sheet
#'
#' The sample sheet must be in TSV format. The first line must be
#' the header, containing column names and optionally starting with
#' a '#' (without any spaces following it).
#'
#' @param samples_sheet String. Path to the TSV file
#' @param column String. Column containing the sample IDs
#' @return Vector of strings. IDs of the samples.
#' @examples
#'   samples <- get_sample_IDs("~/project/samples.tsv", "PatientID")
#'
get_sample_IDs <- function(samples_sheet, column) {
  header <- read.table(samples_sheet, sep="\t", header=FALSE,  nrows=1,
                       stringsAsFactors=FALSE, comment.char="")
  header <- gsub("^#", "", unlist(header))
  samples_metadata <- read.table(samples_sheet, sep="\t", col.names = header,
                                 header=FALSE, skip=1, comment.char="",
                                 colClasses = setNames("character", column))
  print(samples_metadata)
  samples <- subset(samples_metadata, select = c(column))
  samples <- samples[[1]]
  print(paste("Samples:", paste(samples, collapse=", ")))
  samples
}

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
  so$mitoRatio <- PercentageFeatureSet(object = so, pattern = "^MT-")
  so$mitoRatio <- so@meta.data$mitoRatio / 100
  so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
  so$cells <- rownames(so@meta.data)
  so
}

DEFAULT_COUNTS_FILENAME <- "filtered_feature_bc_matrix.h5"

#' Read Counts using Read10X_h5 and call CreateSeuratObject
#'
#' Streamlines reading from multiple directories under a common path
#' containing the output of "Cellranger count" for different samples
#' and adds further QC metrics to the seurat object
#'
#' @param sampleID String. Sample ID of the specific sample
#' @param data_dir_template String. Path to the sample directory, where the
#'                          sample ID is given as '{sampleID}'
#' @param counts_filename String. The filename containing the counts.
#' @param min_cells Integer. The min.cells parameter of CreateSeuratObject
#' @param min_features Integer. The min.features parameter of CreateSeuratObject
#'
#' @return SeuratObject from the given count matrix.
counts_to_seurat <- function(sampleID, data_dir_template,
                          counts_filename = DEFAULT_COUNTS_FILENAME,
                          min_cells = 0, min_features = 0) {
  data_dir <- glue(data_dir_template)
  counts_fullpath <- paste0(data_dir, counts_filename)
  data <- Read10X_h5(filename = counts_fullpath)
  so <- CreateSeuratObject(counts = data,
                           project = sampleID,
                           min.cells = min_cells,
                           min.features = min_features)
  rm(data)
  so <- prepare_QC(so)
  so
}

#' Add QC information from scater library to Seurat object
#'
#' @param so Seurat object
#'
#' @return Seurat object with additional QC metrics
scater_qc <- function(so) {
  library(scater)
  sce <- SingleCellExperiment(assays = list(counts = so@assays$RNA@counts))
  is_mito <- grep("^MT-", rownames(sce))
  qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is_mito))
  filt_ad <- data.frame(
    qc.lib.low = isOutlier(qcstats$sum, log = TRUE, type = "lower"),
    qc.lib.high = isOutlier(qcstats$sum,
                            log = FALSE, type = "higher", nmads = 4),
    qc.nexpr.low = isOutlier(qcstats$detected, log = TRUE, type = "lower"),
    qc.nexpr.high = isOutlier(qcstats$detected,
                              log = FALSE, type = "higher", nmads = 4),
    qc.mito.high = isOutlier(qcstats$subsets_Mito_percent,
                             log = FALSE, type = "higher")
  )
  qcstats <- cbind(qcstats, filt_ad)
  summary(qcstats)
  names(qcstats) <- paste0("scater_", names(qcstats))
  so@meta.data <- cbind(so@meta.data, qcstats)
  so
}

#' Add QC information from robustbase library to Seurat object
#'
#' @param so Seurat object
#'
#' @return Seurat object with additional QC metrics
robustbase_qc <- function(so) {
  library(robustbase)
  set.seed(1234)
  stats <- cbind(log10(so$nCount_RNA), log10(so$nFeature_RNA), so$mitoRatio)
  outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  multi.outlier <- isOutlier(outlying, type = "higher")
  summary(multi.outlier)
  so@meta.data <- cbind(so@meta.data, multi.outlier)
  so
}

#' Identify cell cycle phases
#'
#' derived from a function of Sebastien Mella
#'
id_cc_phase <- function(so) {
  s_genes = cc.genes.updated.2019$s.genes
  g2m_genes = cc.genes.updated.2019$g2m.genes
  s_genes_table <- as.data.frame(table(s_genes %in% rownames(so)))
  g2m_genes_table <- as.data.frame(table(g2m_genes %in% rownames(so)))
  so <- NormalizeData(so, assay = "RNA")
  so <- CellCycleScoring(so, s.features = s_genes[s_genes %in% rownames(so)],
                         g2m.features = g2m_genes[g2m_genes %in% rownames(so)])
  ccp_df <- as.data.frame.array(table(so$Phase))
  ccp_df$phase <- rownames(ccp_df)
  so <- FindVariableFeatures(so, selection.method = "vst", array = "RNA",
                             nfeatures = 2000, verbose = FALSE)
  so <- ScaleData(so, assay = "RNA", verbose = FALSE)
  so <- RunPCA(so, assay = "RNA", reduction.name = "pcaRNA",
               reduction.key = "pcaRNA_", verbose = FALSE)
  ccg_colors <-
    RColorBrewer::brewer.pal(length(levels(so$Phase)), name = "Set1")
  plot(DimPlot(so, reduction = "pcaRNA", group.by= "Phase", cols = ccg_colors))
  plot(DimPlot(so, reduction = "pcaRNA", group.by= "Phase", split.by = "Phase",
               cols = ccg_colors))
  so
}

#' Create Seurat objects for multiple samples
#'
#' the sample IDs are added in the meta.data as 'sample' column
#'
#' @param sampleIDs list of sampleIDs
#' @param input_dir_template String. Path to the sample directory, where the
#'                          sample ID is given as '{sampleID}'
#' @param counts_filename String. The filename containing the counts.
#'
#' @return list of Seurat objects for each of the samples
#'
multi_counts_to_seurat <- function(sampleIDs, input_dir_template,
                                   counts_filename = DEFAULT_COUNTS_FILENAME) {
  so_list <- list()
  for (sampleID in sampleIDs) {
    so_list[[sampleID]] <-
      counts_to_seurat(sampleID, input_dir_template, counts_filename)
    so_list[[sampleID]]$sample <- sampleID
    so_list[[sampleID]] <- scater_qc(so_list[[sampleID]])
    so_list[[sampleID]] <- robustbase_qc(so_list[[sampleID]])
    so_list[[sampleID]] <- id_cc_phase(so_list[[sampleID]])
    print(paste0("Created Seurat object for sample: '", sampleID,
                 "' (n_cells: ", length(so_list[[sampleID]]$cells),
                 ")"))
  }
  so_list
}

#' Sample IDs ordered by number of cells
#'
#' Get the order of sample IDs according to the number of cells
#' in the Seurat objects
#'
#' @param so_list List of Seurat objects
#'
#' @return Vector of strings. Sample IDs in the order of the number of cells
#'
get_n_cells_order <- function(so_list) {
  so_list <- so_list[order(sapply(so_list, function(x) -length(x$cells)))]
  result <- names(so_list)
  print(paste("Sample IDs by decr. n_cells:",
              paste(result, collapse=", ")))
  result
}

#' Merge multiple Seurat objects
#'
#' Merges all Seurat objects in the given list into a single object. using the
#' merge function of SeuratObject. Optionally the sample IDs are reordered. By
#' default the memory of the input list of objects is freed.
#'
#' @param so_list List of Seurat objects
#' @param order Vector of strings. Order of the samples in the merged object.
#'              Defaults to the order of the samples in the list.
#' @param freemem Boolean. Whether to free the memory of the input objects
#'
#' @return A Seurat object.
#'     The cell IDs get the sample IDs prefixed to their barcode.
#'     The sample IDs are transformed into factors to keep their order in plots.
#'
merge_seurat_list <- function(so_list, order=NULL, freemem=TRUE) {
  if (is.null(order)) {
    sampleIDs <- names(so_list)
  }
  else {
    so_list <- so_list[order]
    sampleIDs <- order
  }
  so <- merge(x = so_list[[1]],
              y = so_list[2:length(so_list)],
              add.cell.id = sampleIDs)
  so@meta.data$sample <- factor(so@meta.data$sample, levels = sampleIDs)
  if (freemem) {
    rm(so_list)
  }
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
  ggplot(so@meta.data, aes_string(x=measure)) +
    geom_histogram(bins=100) +
  ggtitle(paste(measure, "Distribution")) +
  xlab(measure) + ylab("Frequency")
}


#' Plot number of cells in different samples
#'
#' @param so Seurat object
#'
histogram_n_cells <- function(so) {
  so@meta.data %>%
    	ggplot(aes(x=sample, fill=sample)) +
    	geom_bar() +
    	theme_classic() +
    	theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    	theme(plot.title = element_text(hjust=0.5, face="bold")) +
    	ggtitle("Number of cells")
}

#' Plot number of UMIs per cell in different samples
#'
#' @param so Seurat object
#'
density_plot_n_umis <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot(aes(color=sample, x=nCount_RNA, fill=sample)) +
    	geom_density(alpha = 0.2) +
    	scale_x_log10() +
    	theme_classic() +
    	ylab("cell density") +
    	geom_vline(xintercept = 100) +
    	geom_vline(xintercept = 500) +
    	geom_vline(xintercept = 1000) +
      facet_wrap(~sample, nrow=plot_nrows) +
    	ggtitle("Number of UMIs per cell")
}

#' Boxplot of the number of UMIs per cell
#'
#' @param so Seurat object
#'
boxplot_n_umis <- function(so) {
  so@meta.data %>%
    	ggplot(aes(x=sample, y=log10(nCount_RNA), fill=sample)) +
    	geom_boxplot() +
    	theme_classic() +
    	theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    	theme(plot.title = element_text(hjust=0.5, face="bold")) +
    	ggtitle("Number of UMIs per cell")
}

#' Plot number of genes per cell in different samples
#'
#' @param so Seurat object
#'
density_plot_n_genes <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot(aes(color=sample, x=nFeature_RNA, fill=sample)) +
    	geom_density(alpha = 0.2) +
    	theme_classic() +
    	scale_x_log10() +
    	geom_vline(xintercept = 300) +
      facet_wrap(~sample, nrow=plot_nrows) +
    	ggtitle("Number of genes per cell")
}


#' Boxplot of the number of genes per cell
#'
#' @param so Seurat object
#'
boxplot_n_genes <- function(so) {
  so@meta.data %>%
    	ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) +
    	geom_boxplot() +
    	theme_classic() +
    	theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    	theme(plot.title = element_text(hjust=0.5, face="bold")) +
    	ggtitle("Number of genes per cell")
}

dotplot_n_umis_genes_mito <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) +
    	geom_point() +
  	  scale_colour_gradient(low = "gray90", high = "black") +
    	stat_smooth(method=lm) +
    	scale_x_log10() +
    	scale_y_log10() +
    	theme_classic() +
    	geom_vline(xintercept = 500) +
    	geom_hline(yintercept = 250) +
    	facet_wrap(~sample, nrow=plot_nrows) +
      ggtitle("N. UMIs vs N. genes and mit.genes ratio")
}

density_plot_mito_ratio <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
    	geom_density(alpha = 0.2) +
    	scale_x_log10() +
    	theme_classic() +
    	geom_vline(xintercept = 0.2) +
      facet_wrap(~sample, nrow=plot_nrows) +
      ggtitle("Mitochondrial gene ratio per cell")
}

density_plot_complexity <- function(so, plot_nrows=0) {
  so@meta.data %>%
    	ggplot(aes(x=log10GenesPerUMI, color=sample, fill=sample)) +
    	geom_density(alpha = 0.2) +
    	theme_classic() +
    	geom_vline(xintercept = 0.8) +
      facet_wrap(~sample, nrow=plot_nrows) +
      ggtitle("Transcriptional complexity")
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

#' Get the path to a step directory
#'
#' @param steps_dir String. Path to the "steps" directory.
#' @param step_n String. Prefix of the step.
#'
#' @return String. Path to the step directory.
get_step_dir <- function(steps_dir, step_n) {
  all_files <- list.files(path = steps_dir)
  matching_files <- all_files[grepl(paste0("^", step_n), all_files)]
  if (length(matching_files) == 0) {
    stop(paste0("No entries matching '",
                step_n, "' in '", steps_dir, "'"))
  }
  if (length(matching_files) > 1) {
    stop(paste0("Multiple entries matching '",
                step_n, "' in '", steps_dir, "'"))
  }
  paste0(steps_dir, "/", matching_files[1], "/")
}

#' Plot the highest expressed genes
#'
#' based on a similar function by Sebastien Mella
#'
#' @param so Seurat object
#' @param n_genes Integer. Number of genes to plot
#' @param sample String. Keep only specified sample. Default: keep all.
#'
#' @return A ggplot object
h_expr_genes_plot <- function(so, n_genes = 10, sample = NULL) {
  library(reshape2)
  library(RColorBrewer)
  library(ggsci)
  library(matrixStats)
  title <- "Highest expressed genes"
  count_matrix <- as.matrix(so@assays$RNA@counts)
  if (!is.null(sample)) {
    count_matrix <- count_matrix[, grepl(sample, colnames(count_matrix))]
    title <- paste0(title, " (Sample ", sample, ")")
  }
  values <- head(order(rowSums2(count_matrix), decreasing = TRUE), n_genes)
  sub_mat <- as.data.frame(count_matrix[values, ])
  sub_mat$feature <- factor(rownames(sub_mat), levels = rev(rownames(sub_mat)))
  sub_mat_long_fmt <- melt(sub_mat, id.vars = "feature")
  plot_colors <- colorRampPalette(
                            rev(pal_futurama("planetexpress")(12)))(n_genes)
  ggplot(sub_mat_long_fmt, aes(x = value, y = feature))+
    geom_boxplot(aes(fill = feature), alpha = .75) +
    scale_fill_manual(values = plot_colors) + theme_light() +
    guides(fill="none") + xlab("Count") + ylab("Gene") + ggtitle(title)
}
