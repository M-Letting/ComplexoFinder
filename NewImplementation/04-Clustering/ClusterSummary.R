# Libraries
library(data.table)


#' Build a cluster summary table from a VSClust membership matrix
#'
#' Splits the underscore-delimited row names of a VSClust membership matrix
#' into separate metadata columns and column-binds them with the membership
#' values, returning one row per peptide feature.
#'
#' Row names are expected to follow the format produced by
#' `prepare_for_vsclust()`: identifier fields pasted together with `"_"`, e.g.
#' `"SMARCA5_P12345_P12345 [10-20]_NA_Phospho"`. A numeric suffix that makes an
#' identifier unique (`".1"`) is removed from the datatype field.
#'
#' @param membership_matrix Numeric matrix. VSClust membership matrix with
#'   peptide identifiers as row names and one column per cluster (e.g.
#'   `vsclust_result$ClustOut$Bestcl$membership`).
#' @param id_field_names Character vector. Names to give the identifier
#'   fields, in the order they appear in the row names. Must not exceed the
#'   number of `"_"`-delimited fields in the row names, and must include
#'   `"Datatype"`.
#'
#' @return A `data.table` with `length(id_field_names)` metadata columns
#'   followed by one membership column per cluster. Number of rows equals
#'   `nrow(membership_matrix)`.
#'
#' @export
make_cluster_summary <- function(
  membership_matrix,
  id_field_names = c("Gene name", "Accession", "Position", "Modification", "Datatype")
) {
  id_fields <- data.table::tstrsplit(
    rownames(membership_matrix),
    "_",
    fixed = TRUE,
    keep = seq_along(id_field_names)
  )
  metadata <- setNames(as.data.table(id_fields), id_field_names)
  metadata[, Datatype := sub("\\..*$", "", Datatype)]

  cbind(metadata, as.data.table(membership_matrix))
}
