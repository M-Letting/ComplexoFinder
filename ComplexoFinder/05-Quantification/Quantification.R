library(data.table)

#' Quantify complexoform-group abundance from peptide intensities and memberships
#'
#' Computes abundance profiles for each complexoform group (membership column)
#' by combining peptide-level intensity values with peptide membership weights.
#'
#' @param peptide_data A `data.table` (or coercible) with peptide identifiers and
#'   intensity columns.
#' @param membership_matrix A numeric matrix/data.frame with peptides in rows and
#'   complexoform groups in columns.
#' @param intensity_cols Character vector of intensity columns. If `NULL`, columns
#'   matching `"^(H9|IMR90)_B[0-9]+_D[0-9]+$"` are used.
#' @param id_cols Character vector used to build peptide identifiers when
#'   `Identifier` is not present.
#' @param aggregation Either `"weighted_mean"` (default) or `"weighted_sum"`.
#' @param membership_threshold Numeric threshold; memberships below this are set
#'   to zero before quantification.
#' @param normalize_membership Logical. If `TRUE`, memberships are row-normalized
#'   to sum to 1 after thresholding.
#'
#' @return A list with:
#'   - `abundance_wide`: one row per complexoform group, one column per sample
#'   - `group_stats`: number of contributing peptides and total membership weight
#'
#' @export
quantify_complexoform_abundance <- function(
  peptide_data,
  membership_matrix,
  intensity_cols = NULL,
  id_cols = c("Gene name", "Peptide"),
  aggregation = c("weighted_mean", "weighted_sum"),
  membership_threshold = 0,
  normalize_membership = FALSE,
  cond_regex = NULL
) {
  # Resolve aggregation option early so downstream branches are explicit.
  aggregation <- match.arg(aggregation)

  # Coerce to data.table to keep subsetting and joins consistent.
  if (!is.data.table(peptide_data)) {
    peptide_data <- as.data.table(peptide_data)
  }

  # Auto-detect intensity columns when not provided by the caller.
  if (is.null(intensity_cols)) {
    intensity_cols <- grep(
      "^(H9|IMR90)_B[0-9]+_D[0-9]+$",
      names(peptide_data),
      value = TRUE
    )
  }

  if (length(intensity_cols) == 0L) {
    stop("No intensity columns found. Provide intensity_cols explicitly.")
  }

  missing_intensity <- setdiff(intensity_cols, names(peptide_data))
  if (length(missing_intensity) > 0L) {
    stop(
      "Missing intensity columns: ",
      paste(missing_intensity, collapse = ", ")
    )
  }

  if (is.null(rownames(membership_matrix))) {
    stop("membership_matrix must have rownames matching peptide identifiers.")
  }

  membership_matrix <- as.matrix(membership_matrix)
  storage.mode(membership_matrix) <- "numeric"

  # Build peptide identifiers either from an existing Identifier field
  # or from user-specified identifier columns.
  if ("Identifier" %in% names(peptide_data)) {
    peptide_ids <- peptide_data$Identifier
  } else {
    missing_id_cols <- setdiff(id_cols, names(peptide_data))
    if (length(missing_id_cols) > 0L) {
      stop(
        "Missing identifier columns: ",
        paste(missing_id_cols, collapse = ", ")
      )
    }
    peptide_ids <- do.call(paste, c(peptide_data[, ..id_cols], sep = "_"))
  }

  # VSClust may disambiguate duplicated feature identifiers using make.unique()
  # (e.g., appending .1, .2). Mirror that behavior here so alignment is robust.
  n_dup_ids <- anyDuplicated(peptide_ids)
  if (n_dup_ids > 0L) {
    peptide_ids <- make.unique(as.character(peptide_ids))
  }

  common_ids <- intersect(peptide_ids, rownames(membership_matrix))
  if (length(common_ids) == 0L) {
    stop(
      "No overlapping peptide identifiers between peptide_data and memberships."
    )
  }

  # Align peptide intensity rows and membership rows to the same ID order.
  peptide_aligned <- peptide_data[match(common_ids, peptide_ids)]
  membership_aligned <- membership_matrix[common_ids, , drop = FALSE]

  if (!is.numeric(membership_threshold) || length(membership_threshold) != 1L) {
    stop("membership_threshold must be a single numeric value.")
  }

  membership_aligned[membership_aligned < membership_threshold] <- 0

  # Optionally normalize each peptide's memberships to sum to 1.
  if (normalize_membership) {
    row_sums <- rowSums(membership_aligned, na.rm = TRUE)
    valid_rows <- row_sums > 0
    membership_aligned[valid_rows, ] <- membership_aligned[
      valid_rows,
      ,
      drop = FALSE
    ] /
      row_sums[valid_rows]
  }

  intensity_mat <- as.matrix(peptide_aligned[, ..intensity_cols])
  storage.mode(intensity_mat) <- "numeric"

  # Pre-allocate output matrix: rows = groups, cols = samples.
  n_groups <- ncol(membership_aligned)
  n_samples <- ncol(intensity_mat)
  abundance_mat <- matrix(
    NA_real_,
    nrow = n_groups,
    ncol = n_samples,
    dimnames = list(colnames(membership_aligned), colnames(intensity_mat))
  )

  for (group_idx in seq_len(n_groups)) {
    w <- membership_aligned[, group_idx]

    # Weighted sum: preserves total signal contribution per sample.
    if (aggregation == "weighted_sum") {
      abundance_mat[group_idx, ] <- colSums(intensity_mat * w, na.rm = TRUE)
    } else {
      # Weighted mean: averages using only non-missing peptide intensities.
      numerator <- colSums(intensity_mat * w, na.rm = TRUE)
      denominator <- colSums((!is.na(intensity_mat)) * w, na.rm = TRUE)
      vals <- numerator / denominator
      vals[denominator == 0] <- NA_real_
      abundance_mat[group_idx, ] <- vals
    }
  }

  group_labels <- rownames(abundance_mat)
  if (is.null(group_labels)) {
    group_labels <- paste0("Group", seq_len(nrow(abundance_mat)))
  }

  # Optionally collapse replicate columns to per-condition means.
  # cond_regex must have one capture group that extracts the condition label.
  if (!is.null(cond_regex)) {
    sample_names <- colnames(abundance_mat)
    cond_labels <- sub(cond_regex, "\\1", sample_names, perl = TRUE)
    unique_conds <- unique(cond_labels)

    cond_mat <- matrix(
      NA_real_,
      nrow = nrow(abundance_mat),
      ncol = length(unique_conds),
      dimnames = list(group_labels, unique_conds)
    )
    for (cond in unique_conds) {
      cols <- which(cond_labels == cond)
      cond_mat[, cond] <- rowMeans(
        abundance_mat[, cols, drop = FALSE],
        na.rm = TRUE
      )
    }
    cond_mat[is.nan(cond_mat)] <- NA_real_

    abundance_wide <- data.table(
      ComplexoformGroup = group_labels,
      as.data.table(cond_mat)
    )
  } else {
    abundance_wide <- data.table(
      ComplexoformGroup = group_labels,
      as.data.table(abundance_mat)
    )
  }

  # Group-level diagnostics for interpretation and QC.
  group_stats <- data.table(
    ComplexoformGroup = colnames(membership_aligned),
    NPeptides = colSums(membership_aligned > 0, na.rm = TRUE),
    MembershipWeightSum = colSums(membership_aligned, na.rm = TRUE)
  )

  list(
    abundance_wide = abundance_wide,
    group_stats = group_stats
  )
}


#' Quantify complexoform abundance from a VSClust result object
#'
#' Convenience wrapper around `quantify_complexoform_abundance()` using
#' `vsclust_result$original_data` and
#' `vsclust_result$ClustOut$Bestcl$membership`.
#'
#' @inheritParams quantify_complexoform_abundance
#' @param vsclust_result Result object from `vsclust_on_complex()`.
#'
#' @return Same as `quantify_complexoform_abundance()`.
#'
#' @export
quantify_complexoform_abundance_from_vsclust <- function(
  vsclust_result,
  intensity_cols = NULL,
  id_cols = c("Gene name", "Peptide"),
  aggregation = c("weighted_mean", "weighted_sum"),
  membership_threshold = 0,
  normalize_membership = FALSE,
  cond_regex = NULL
) {
  # Extract required objects from the VSClust pipeline result.
  membership_matrix <- vsclust_result$ClustOut$Bestcl$membership
  peptide_data <- vsclust_result$original_data

  if (is.null(membership_matrix)) {
    stop("vsclust_result$ClustOut$Bestcl$membership is missing.")
  }
  if (is.null(peptide_data)) {
    stop("vsclust_result$original_data is missing.")
  }

  # Delegate to the core quantification routine.
  quantify_complexoform_abundance(
    peptide_data = peptide_data,
    membership_matrix = membership_matrix,
    intensity_cols = intensity_cols,
    id_cols = id_cols,
    aggregation = aggregation,
    membership_threshold = membership_threshold,
    normalize_membership = normalize_membership,
    cond_regex = cond_regex
  )
}
