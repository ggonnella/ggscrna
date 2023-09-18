
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
