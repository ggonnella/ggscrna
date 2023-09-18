
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
#' samples <- GetSampleIDs("~/project/samples.tsv", "PatientID")
#'
GetSampleIDs <- function(samples_sheet, column) {
  header <- read.table(samples_sheet, sep="\t", header=FALSE,  nrows=1,
                       stringsAsFactors=FALSE, comment.char="")
  header <- gsub("^#", "", unlist(header))
  samples_metadata <- read.table(samples_sheet, sep="\t", col.names = header,
                                 header=FALSE, skip=1, comment.char="",
                                 colClasses = setNames("character", column))
  print(samples_metadata)
  samples <- subset(samples_metadata, select = c(column))
  samples[[1]]
}

#' Read Counts using Read10X_h5 and call CreateSeuratObject
#'
#' Streamlines reading from multiple directories under a common path
#' containing the output of "Cellranger count" for different samples
#'
#' @param sampleID String. Sample ID of the specific sample
#' @param data_dir_template String. Path to the sample directory, where the
#'                          sample ID is given as '{sampleID}'
#' @param counts_filename String. The filename containing the counts.
#' @param min_cells Integer. The min.cells parameter of CreateSeuratObject
#' @param min_features Integer. The min.features parameter of CreateSeuratObject
#'
#' @returns SeuratObject from the given count matrix.
Counts2Seurat <- function(sampleID, data_dir_template,
                          counts_filename = "filtered_feature_bc_matrix.h5",
                          min_cells = 0, min_features = 0) {
  data_dir <- glue(data_dir_template)
  counts_fullpath <- paste0(data_dir, counts_filename)
  print(counts_fullpath)
  data <- Read10X_h5(filename = counts_fullpath)
  so <- CreateSeuratObject(counts = data,
                           project = sampleID,
                           min.cells = min_cells,
                           min.features = min_features)
  rm(data)
  so
}

#' Compute additional metrics on the Seurat object
#'
#' The metrics are those explained in:
#' https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
#'
#' ## Notes
#' - mt genes computation by MT- prefix assumes this is a human dataset,
#'   hence the name
#'
#' @param so SeuratObject for which to compute the metrics
#'
#' @return The SeuratObject, with additional metrics in ``@meta.data``
ComputeQCMetricsHuman <- function(so) {
  so$mitoRatio <- PercentageFeatureSet(object = so, pattern = "^MT-")
  so$mitoRatio <- so@meta.data$mitoRatio / 100
  so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
  so@meta.data <- so@meta.data %>%
    dplyr::rename(nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  so$cells <- rownames(so@meta.data)
  so
}

#' Create Seurat objects for multiple samples
#'
#' Creates multiple Seurat objects for a list of samples
#' and adds further QC metrics to them
#'
#' @param sampleIDs list of sampleIDs
#' @param data_dir_template String. Path to the sample directory, where the
#'                          sample ID is given as '{sampleID}'
#' @param counts_filename String. The filename containing the counts.
#'
#' @returns list of Seurat objects for each of the samples
#'
SeuratReadMulti <- function(sampleIDs, data_dir_template,
                            counts_filename = "filtered_feature_bc_matrix.h5") {
  so_list <- list()
  so <- NULL
  cell_id_pfx <- list()
  print(samples)
  for (sampleID in samples) {
    so <- Counts2Seurat(sampleID, data_dir_template, counts_filename)
    so <- ComputeQCMetricsHuman(so)
    #print(summary(so@meta.data))
    so_list[[sampleID]] <- so
    rm(so)
    so_list[[sampleID]]$sample <- sampleID
  }
  so_list
}

#' Merge multiple Seurat objects
#'
#' Merges all Seurat objects in the given list into a single object.
#' It is a wrapper for the merge function of SeuratObject.
#'
#' @param so_list List of Seurat objects
#' @param sampleIDs List of strings, Sample IDs in the same order as in so_list
#'
#' @returns A Seurat object.
#'     The cell IDs get the sample name prefixed to their barcode.
#'
MergeMultiSeurat <- function(so_list, sampleIDs) {
  merge(x = so_list[[1]],
        y = so_list[2:length(so_list)],
        add.cell.id = sampleIDs)
}
