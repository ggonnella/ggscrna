#'
#' Giorgio's functions for single cell sequencing analysis
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

#'
#' Acknowledgements:
#' some of the code was originally derived from code written by Sebastien Mella
#'

#' Get sample IDs from a sample sheet
#'
#' The sample sheet must be in TSV format. The first line must be
#' the header, containing column names and optionally starting with
#' a '#' (without any spaces following it).
#'
#' @param samples_sheet String. Path to the TSV file
#' @param column String. Column containing the sample IDs
#' @param verbose Boolean. Whether to print the sample sheet and the IDs.
#' @return Vector of strings. IDs of the samples.
#' @examples
#'   samples <- get_sample_IDs("~/project/samples.tsv", "PatientID")
#'
get_sample_IDs <- function(samples_sheet, column, verbose=TRUE) {
  header <- read.table(samples_sheet, sep="\t", header=FALSE,  nrows=1,
                       stringsAsFactors=FALSE, comment.char="")
  header <- gsub("^#", "", unlist(header))
  samples_metadata <- read.table(samples_sheet, sep="\t", col.names = header,
                                 header=FALSE, skip=1, comment.char="",
                                 colClasses = setNames("character", column))
  if (verbose) {
    print(paste("Sample sheet:", samples_sheet))
    print(samples_metadata)
  }
  samples <- subset(samples_metadata, select = c(column))
  samples <- samples[[1]]
  if (verbose) {
    print(paste("Samples:", paste(samples, collapse=", ")))
  }
  samples
}

CELLRANGER_FILTERED_COUNTS <- "filtered_feature_bc_matrix.h5"

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
                          counts_filename = CELLRANGER_FILTERED_COUNTS,
                          min_cells = 0, min_features = 0) {
  data_dir <- glue::glue(data_dir_template)
  counts_fullpath <- paste0(data_dir, counts_filename)
  data <- Seurat::Read10X_h5(filename = counts_fullpath)
  so <- Seurat::CreateSeuratObject(counts = data,
                           project = sampleID,
                           min.cells = min_cells,
                           min.features = min_features)
  rm(data)
  so <- prepare_QC(so)
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
                                   counts_filename =
                                     CELLRANGER_FILTERED_COUNTS) {
  so_list <- list()
  for (sampleID in sampleIDs) {
    so_list[[sampleID]] <-
      counts_to_seurat(sampleID, input_dir_template, counts_filename)
    so_list[[sampleID]]$sample <- sampleID
    so_list[[sampleID]] <- scater_qc(so_list[[sampleID]])
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
              merge.data = TRUE,
              add.cell.id = sampleIDs)
  so@meta.data$sample <- factor(so@meta.data$sample, levels = sampleIDs)
  if (freemem) {
    rm(so_list)
  }
  so
}

