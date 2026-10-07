#' Find Number of Clusters from dCF Labels
#'
#' Counts the number of unique dCF labels in the data, excluding the
#' "dCF-1" label which represents singletons.
#' This function is used to determined number of clusters for VSClust.
#'
#' @param data data.table containing a "dCF" column with cluster labels
#'
#' @return Integer number of clusters (unique dCF labels excluding "dCF-1"
#'
#' @details The function checks for the presence of a "dCF" column and counts
#' the unique values in that column, excluding the "dCF-1" label.
#'
#' @export

find_nclust <- function(data) {
  # Sanity check
  # Accept data.frames and data.tables; coerce data.frames to data.table
  if (!data.table::is.data.table(data) && !is.data.frame(data)) {
    stop("Input data must be a data.table or data.frame.")
  }
  if (!data.table::is.data.table(data) && is.data.frame(data)) {
    data <- data.table::as.data.table(data)
  }
  if (!"dCF" %in% colnames(data)) {
    stop("Input data must contain a 'dCF' column.")
  }

  # Find number of unique dCF values that are not "dCF-1"
  nclust <- length(unique(data$dCF[data$dCF != "dCF-1"]))
  return(nclust)
}

#' Prepare Data for VSClust
#'
#' Formats data for VSClust analysis by creating unique row identifiers
#' and extracting measurement columns.
#'
#' @param data data.table with proteomics data
#' @param n_rep Integer. Number of replicates
#' @param n_cond Integer. Number of conditions/time points
#' @param id_cols Character vector of columns to create row identifiers
#' @param measurement_cols Character vector of measurement column names
#'
#' @return data.table formatted for VSClust with Identifier column
#'
#' @details
#' The function expects measurement_cols to have length n_rep * n_cond.
#' Row identifiers are created by pasting id_cols together with "_".
#'
#' @export
prepare_for_vsclust <- function(
  data,
  n_rep = 3,
  n_cond = 5,
  id_cols = NULL,
  measurement_cols = NULL
) {
  # Sanity check
  if (is.null(id_cols) || is.null(measurement_cols)) {
    stop("id_cols and measurement_cols must be provided")
  }
  if (length(measurement_cols) != n_rep * n_cond) {
    stop("Amount of measurement columns must equal n_rep * n_cond")
  }

  # Ensure data is a data.table
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }

  # Convert id_cols to unique row identifiers
  row_ids <- apply(data[, id_cols, with = FALSE], 1, function(row) {
    paste(row, collapse = "_")
  })
  row_ids <- make.unique(row_ids)

  # Extract measurement data
  measurement_data <- data[, measurement_cols, with = FALSE]

  # Combine row_ids with measurement data
  cbind(Identifier = row_ids, measurement_data)
}


#' Reorder Columns for VSClust
#'
#' Reorders measurement columns from condition-grouped to replicate-grouped
#' format as expected by VSClust, and converts to data.frame with row names.
#'
#' @param data data.table from prepare_for_vsclust
#' @param n_rep Integer. Number of replicates
#' @param n_cond Integer. Number of conditions
#' @param id_col Character. Column name containing identifiers (default: "Identifier")
#' @param reorder Logical. Whether to reorder columns (default: TRUE)
#'
#' @return data.frame with measurement values and row identifiers as row names
#'
#' @details
#' VSClust expects data in replicate-grouped format:
#' Rep1_Cond1, Rep1_Cond2, ..., Rep2_Cond1, Rep2_Cond2, ...
#'
#' If your data is condition-grouped (Cond1_Rep1, Cond1_Rep2, ...),
#' set reorder = TRUE.
#'
#' @export
reorder_for_vsclust <- function(
  data,
  n_rep = 3,
  n_cond = 5,
  id_col = "Identifier",
  reorder = TRUE
) {
  # Sanity check
  if (!(id_col %in% names(data))) {
    stop(paste("id_col", id_col, "not found in data"))
  }

  # Save row identifiers
  row_ids <- data[[id_col]]
  data <- data[, -..id_col]

  # Check dimensions
  if (n_rep * n_cond != ncol(data)) {
    stop("Number of columns in data does not match n_rep * n_cond")
  }

  if (reorder) {
    # Create index matrix (conditions × replicates)
    idx_matrix <- matrix(
      seq_len(n_rep * n_cond),
      nrow = n_cond,
      ncol = n_rep,
      byrow = TRUE
    )

    # Transpose and flatten to get new column order
    new_order <- as.vector(idx_matrix)

    # Return reordered data.table
    data <- data[, ..new_order]
  }

  # Convert to data.frame for vsclust
  data <- as.data.frame(data)

  # Restore row names
  row.names(data) <- row_ids

  data
}

#' Run VSClust on a single complex-level peptide table
#'
#' Prepares complex-specific peptide intensity data, optionally applies peptide-
#' to-cluster restriction constraints, and runs VSClust clustering with a fixed
#' number of clusters.
#'
#' @param data A `data.table` (or `data.frame`) containing peptide metadata and
#'   intensity columns for one complex.
#' @param id_cols Character vector of column names used to build a unique
#'   peptide identifier for VSClust input (for example protein, gene, peptide,
#'   and modification columns).
#' @param intensity_cols Character vector of intensity columns (columns to be
#'   used for clustering).
#' @param n_rep Integer. Number of replicate measurements per condition.
#' @param n_cond Integer or `NULL`. Number of conditions. If `NULL`,
#'   inferred as `length(intensity_cols) / n_rep`.
#' @param n_clusters Integer. Number of clusters (`NClust`) to fit in VSClust.
#' @param grouped_replicates Logical. Controls column ordering before VSClust.
#'   Use `TRUE` when input is condition-grouped and should be reordered to
#'   replicate-grouped format expected by VSClust.
#' @param restriction_matrix Optional logical matrix encoding peptide-to-cluster
#'   restrictions. Rows must correspond to peptide identifiers used in VSClust
#'   input. `TRUE` values indicate restricted/forbidden assignments in the
#'   constrained clustering step.
#' @param sds_amplification_factor Numeric. Reserved tuning parameter for
#'   standard-deviation amplification in custom workflows.
#'
#' @return A named list with:
#' \describe{
#'   \item{`ClustOut`}{VSClust output object from
#'   `runClustWrapper_custom()`.}
#'   \item{`original_data`}{Original input data (post coercion to
#'   `data.table`).}
#' }
#'
#' @details
#' Workflow performed by this function:
#' \enumerate{
#'   \item Validate inputs and infer `n_cond` when missing.
#'   \item Build unique row identifiers from `id_cols` and extract
#'   `intensity_cols`.
#'   \item Reorder measurement columns for VSClust layout when requested.
#'   \item Run pre-clustering statistics via `PrepareForVsClust_custom()`.
#'   \item Run constrained or unconstrained clustering via
#'   `runClustWrapper_custom()` depending on `restriction_matrix`.
#' }
#'
#' If a `restriction_matrix` is supplied, rows are aligned to peptide
#' identifiers. Missing peptides are automatically added with no restrictions,
#' and non-overlapping matrices are ignored with a warning.
#'
#' @export
vsclust_on_complex <- function(
  data,
  id_cols,
  intensity_cols,
  n_rep = 3,
  n_cond = NULL,
  n_clusters = NULL,
  grouped_replicates = TRUE,
  restriction_matrix = NULL,
  sds_amplification_factor = 1
) {
  # Sanity checks
  if (!requireNamespace("vsclust", quietly = TRUE)) {
    stop(
      "Package 'VSClust' is required. Install with: BiocManager::install('VSClust')"
    )
  }
  if (is.null(n_cond)) {
    n_cond <- length(intensity_cols) / n_rep
  }
  if (length(intensity_cols) != n_cond * n_rep) {
    stop("Length of intensity_cols must equal n_cond * n_rep")
  }
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  if (is.null(n_clusters)) {
    stop("n_clusters must be provided.")
  }

  # Check for sufficient data after filtering
  if (nrow(data) < 3) {
    stop(
      "Too few rows remaining after filtering (",
      nrow(data),
      "). Cannot proceed with clustering."
    )
  }

  # Prepare data for vsclust
  message("Preparing data for vsclust...")
  dat1 <- prepare_for_vsclust(
    data,
    n_rep = n_rep,
    n_cond = n_cond,
    id_cols = id_cols,
    measurement_cols = intensity_cols
  )

  dat1 <- reorder_for_vsclust(
    dat1,
    n_rep = n_rep,
    n_cond = n_cond,
    id_col = "Identifier",
    reorder = grouped_replicates
  )

  # Match and reorder restriction matrix to dat1 if provided
  if (!is.null(restriction_matrix)) {
    # Check if rownames match
    common_peptides <- intersect(rownames(dat1), rownames(restriction_matrix))

    if (length(common_peptides) == 0) {
      warning(
        "No matching peptides between data and restriction matrix. Clustering without constraints."
      )
      restriction_matrix <- NULL
    } else if (length(common_peptides) < nrow(dat1)) {
      message(
        "Restriction matrix available for ",
        length(common_peptides),
        " of ",
        nrow(dat1),
        " peptides"
      )

      # Create a new restriction matrix for all peptides in dat1
      # Missing peptides get FALSE (allowed in all clusters)
      full_restriction_matrix <- matrix(
        FALSE,
        nrow = nrow(dat1),
        ncol = ncol(restriction_matrix)
      )
      rownames(full_restriction_matrix) <- rownames(dat1)
      colnames(full_restriction_matrix) <- colnames(restriction_matrix)

      # Fill in restrictions for peptides that exist in restriction_matrix
      full_restriction_matrix[common_peptides, ] <- restriction_matrix[
        common_peptides,
      ]
      restriction_matrix <- full_restriction_matrix
    } else {
      # All peptides match - just reorder to match dat1
      restriction_matrix <- restriction_matrix[rownames(dat1), , drop = FALSE]
    }
  }

  # Set parameters for vsclust
  cores <- parallel::detectCores() - 1

  # Run pre-clustering analysis
  message("Running pre-clustering analysis...")
  suppressWarnings(
    statOut <- PrepareForVsClust_custom(
      dat1,
      n_rep,
      n_cond,
      isPaired = FALSE,
      TRUE
    )
  )
  dat <- statOut$dat

  Sds <- dat[, ncol(dat)]
  cat(paste(
    "Features:",
    nrow(dat),
    "\nMissing values:",
    sum(is.na(dat)),
    "\nMedian standard deviations:",
    round(median(Sds, na.rm = TRUE), digits = 3),
    "\n"
  ))

  # Run clustering with constraints
  if (!is.null(restriction_matrix)) {
    message(
      "Running clustering with ",
      sum(restriction_matrix),
      " restrictions..."
    )
  } else {
    message("Running clustering without constraints...")
  }

  ClustOut <- runClustWrapper_custom(
    dat,
    NClust = n_clusters,
    constraints = restriction_matrix,
    proteins = NULL,
    VSClust = TRUE,
    scaling = "standardize",
    cores = cores
  )

  # Return clustering results and original data (before vsclust processing)
  return(list(
    ClustOut = ClustOut,
    original_data = data
  ))
}
