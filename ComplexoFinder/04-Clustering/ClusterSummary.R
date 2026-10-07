# Libraries
library(data.table)


#' Build a cluster summary table from a VSClust membership matrix
#'
#' Parses the underscore-delimited row names of a VSClust membership matrix
#' into separate metadata columns and column-binds them with the membership
#' values, returning one row per peptide feature.
#'
#' Row names are expected to follow the format produced by
#' `prepare_for_vsclust()`: fields from `id_cols` pasted together with `"_"`.
#' For example, a row name `"SMARCA5_P12345_[10-20]_Phospho_PTM"` with the
#' default `id_cols` would be parsed into five separate columns.
#'
#' @param membership_matrix Numeric matrix. VSClust membership matrix with
#'   peptide identifiers as row names and one column per cluster (e.g.
#'   `vsclust_result$ClustOut$Bestcl$membership`).
#' @param id_cols Character vector. Names to assign to the parsed identifier
#'   fields, in the order they appear in the row names. Must not exceed the
#'   number of `"_"`-delimited tokens in the row names.
#'
#' @return A `data.table` with `length(id_cols)` metadata columns followed by
#'   one membership column per cluster. Number of rows equals
#'   `nrow(membership_matrix)`.
#'
#' @export
make_cluster_summary <- function(
  membership_matrix,
  id_cols = c("Gene name", "Accession", "Position", "Modification", "Datatype")
) {
  parts <- data.table::tstrsplit(
    rownames(membership_matrix),
    "_",
    fixed = TRUE,
    keep = seq_along(id_cols)
  )
  metadata <- setNames(as.data.table(parts), id_cols)
  metadata[, Datatype := sub("\\..*$", "", Datatype)]

  cbind(metadata, as.data.table(membership_matrix))
}
