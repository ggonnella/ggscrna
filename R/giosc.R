#'
#' Functions for reading/merging/splitting a Seurat object.
#'

# Licence: CC-BY-SA
# (c) Giorgio Gonnella, 2023

#'
#' Acknowledgements:
#' some of the code was originally derived from code written by Sebastien Mella
#'

CELLRANGER_FILTERED_COUNTS <- "filtered_feature_bc_matrix.h5"
CELLRANGER_RAW_COUNTS <- "raw_feature_bc_matrix.h5"

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

#' Check that a sample exists and extract it
#'
#' @param so          Seurat object
#' @param samplename  String. Sample name
#' 
#' @return           Seurat object with only the cells of the specified sample
so_extract_sample <- function(so, sample) {
  if (!sample %in% so$sample)
    stop(paste0("Sample ", sample, " not found in Seurat object"))
  subset(so, sample == sample)
}

#' Read multiple Seurat Objects from RDS for each sample in a samples_sheet
#' 
#' @param samples_sheet        String. Path to the samples sheet
#' @param samples_sheet_idcol  String. Column name of the sample ID in the sheet
#' @param file_pfx             String. Prefix of the RDS files
#' @param file_sfx             String. Suffix of the RDS files
#' @param verbose              (optional) Boolean. Whether to print the sample IDs
#'                             Defaults to FALSE
#' 
#' @return                     List of Seurat objects
#' 
read_samples_rds <- function(samples_sheet, samples_sheet_idcol, file_pfx, file_sfx,
                             verbose=FALSE) {
  sampleIDs <- get_sample_IDs(samples_sheet, column=samples_sheet_idcol)
  so_list <- list()
  for (sampleID in sampleIDs) {
    rds_file <- paste0(file_pfx, sampleID, file_sfx)
    so_list[[sampleID]] <- readRDS(rds_file)
    if (verbose)
      print(paste0("Read Seurat object from file '", rds_file, "'for sample: '", sampleID, "'"))
  }
  so_list
}

#' Create a Seurat object for testing
#' 
#' @param num_genes             Integer. Number of genes.
#' @param num_cells_per_sample  Integer. Number of cells per sample.
#' @param num_samples           Integer. Number of samples.
#' @param low_expr_part         Numeric. Proportion of genes with low expression.
#'                                       (between 0 and 1)
#' @param low_avg               Numeric. Average of low expression genes.
#'                                       (used as lambda of Poisson distribution)
#' @param high_avg              Numeric. Average of high expression genes.
#'                                       (used as lambda of Poisson distribution)
#' 
#' @return Seurat object.
#'
create_test_so <- function(num_genes, num_cells_per_sample, num_samples,
                           low_expr_part = 0.8, low_avg = 2, high_avg = 20) {
  total_cells <- num_cells_per_sample * num_samples

  n_low_expr <- round(num_genes * low_expr_part)
  counts_low <- matrix(rpois(n_low_expr * total_cells, low_avg), 
                      nrow = n_low_expr, ncol = total_cells)

  n_high_expr <- num_genes - n_low_expr
  counts_high <- matrix(rpois(n_high_expr * total_cells, high_avg), 
                       nrow = n_high_expr, ncol = total_cells)

  # join and then random shuffle to mix high and low:
  counts <- rbind(counts_low, counts_high)
  counts = counts[sample(nrow(counts)), ]

  genes <- paste0("Gene", seq_len(num_genes))
  rownames(counts) <- genes

  samples <- rep(paste0("Sample", seq_len(num_samples)), each = num_cells_per_sample)
  cell_names <- paste0("Cell", seq_len(total_cells), '_', samples)
  colnames(counts) <- cell_names

  so <- CreateSeuratObject(counts = counts, project = "TestData")
  so$sample <- samples

  so
}

#' Consolidate multiple features of a Seurat object
#' 
#' @param so               Seurat object
#' @param features_to_sum  Vector of strings. The names of the features to sum
#' @param new_feature_name String. The name of the new feature
#' @param assay            String. Name of the assay to use.
#' @param layer            String. Name of the layer to use.
#' 
#' @return                 Seurat object with the specified features summed
#'
#' Note: the consolidation is done by creating a new object to which the counts and metadata
#'       are copied. Only a single assay and layer is copied. The rest is not included in the 
#'       constructed so 
#'
consolidate_features <- function(so, features_to_sum, new_feature_name, assay = "RNA", layer = "counts") {
  matrix <- GetAssayData(object = so, layer = layer, assay = assay)
  new_matrix <- consolidate_matrix_rows(matrix, features_to_sum,
                                        new_feature_name, keep_length = TRUE)
  metadata <- as.data.frame(so@metadata)
  result <- CreateSeuratObject(counts = new_matrix, assay = assay, meta.data = metadata, project = so$project)

  result
}