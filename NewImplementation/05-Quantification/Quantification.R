library(data.table)

#' Quantify complexoform-group abundance from peptide intensities and memberships
#'
#' Computes an abundance profile for each complexoform group (membership
#' column) as the membership-weighted mean (or sum) of the peptide
#' intensities.
#'
#' @param peptide_data A `data.table` (or coercible) with peptide identifiers and
#'   intensity columns.
#' @param membership_matrix A numeric matrix/data.frame with peptides in rows and
#'   complexoform groups in columns.
#' @param intensity_columns Character vector of intensity columns. If `NULL`,
#'   columns matching `"^(H9|IMR90)_B[0-9]+_D[0-9]+$"` are used.
#' @param id_columns Character vector used to build peptide identifiers when
#'   `Identifier` is not present.
#' @param aggregation Either `"weighted_mean"` or `"weighted_sum"`.
#' @param membership_threshold Numeric threshold; memberships below this are set
#'   to zero before quantification, so that only peptides with at least this
#'   membership in a group contribute to it. With 0 every peptide contributes
#'   to every group.
#' @param normalize_membership Logical. If `TRUE`, memberships are row-normalized
#'   to sum to 1 after thresholding.
#' @param center_peptides Logical. If `TRUE`, each peptide's own mean level is
#'   removed before averaging and the group's mean level is added back, so a
#'   group's abundance does not shift when a peptide with a high or low level
#'   is missing from a sample. Without missing values the result is the same
#'   as the plain weighted mean. Only for `"weighted_mean"`.
#' @param condition_regex Optional regex with one capture group for the
#'   condition label; if given, replicate columns are averaged per condition.
#'
#' @return A list with:
#'   - `abundance_wide`: one row per complexoform group, one column per sample
#'     (or per condition, with `condition_regex`)
#'   - `group_stats`: number of contributing peptides and total membership weight
#'
#' @export
quantify_complexoform_abundance <- function(
  peptide_data,
  membership_matrix,
  intensity_columns = NULL,
  id_columns = c("Gene name", "Peptide"),
  aggregation = c("weighted_mean", "weighted_sum"),
  membership_threshold = 0.5,
  normalize_membership = FALSE,
  center_peptides = TRUE,
  condition_regex = NULL
) {
  aggregation <- match.arg(aggregation)

  if (!is.data.table(peptide_data)) {
    peptide_data <- as.data.table(peptide_data)
  }

  # Detect the intensity columns when the caller does not give them
  if (is.null(intensity_columns)) {
    intensity_columns <- grep(
      "^(H9|IMR90)_B[0-9]+_D[0-9]+$",
      names(peptide_data),
      value = TRUE
    )
  }

  if (length(intensity_columns) == 0L) {
    stop("No intensity columns found. Provide intensity_columns explicitly.")
  }

  missing_intensity_columns <- setdiff(intensity_columns, names(peptide_data))
  if (length(missing_intensity_columns) > 0L) {
    stop(
      "Missing intensity columns: ",
      paste(missing_intensity_columns, collapse = ", ")
    )
  }

  if (is.null(rownames(membership_matrix))) {
    stop("membership_matrix must have rownames matching peptide identifiers.")
  }

  membership_matrix <- as.matrix(membership_matrix)
  storage.mode(membership_matrix) <- "numeric"

  # Peptide identifiers: an existing Identifier column, or the identifier
  # columns pasted together
  if ("Identifier" %in% names(peptide_data)) {
    peptide_ids <- peptide_data$Identifier
  } else {
    missing_id_columns <- setdiff(id_columns, names(peptide_data))
    if (length(missing_id_columns) > 0L) {
      stop(
        "Missing identifier columns: ",
        paste(missing_id_columns, collapse = ", ")
      )
    }
    peptide_ids <- do.call(paste, c(peptide_data[, ..id_columns], sep = "_"))
  }

  if (anyDuplicated(peptide_ids) > 0L) {
    stop("Peptide identifiers are not unique.")
  }

  shared_ids <- intersect(peptide_ids, rownames(membership_matrix))
  if (length(shared_ids) == 0L) {
    stop(
      "No overlapping peptide identifiers between peptide_data and memberships."
    )
  }

  # Peptide rows and membership rows in the same order
  aligned_peptides <- peptide_data[match(shared_ids, peptide_ids)]
  membership <- membership_matrix[shared_ids, , drop = FALSE]

  if (!is.numeric(membership_threshold) || length(membership_threshold) != 1L) {
    stop("membership_threshold must be a single numeric value.")
  }

  membership[membership < membership_threshold] <- 0

  # Optionally normalize each peptide's memberships to sum to 1
  if (normalize_membership) {
    membership_total <- rowSums(membership, na.rm = TRUE)
    has_membership <- membership_total > 0
    membership[has_membership, ] <- membership[
      has_membership,
      ,
      drop = FALSE
    ] /
      membership_total[has_membership]
  }

  intensities <- as.matrix(aligned_peptides[, ..intensity_columns])
  storage.mode(intensities) <- "numeric"

  # Each peptide's own level (its mean over the samples it was measured in)
  peptide_level <- rowMeans(intensities, na.rm = TRUE)
  if (center_peptides && aggregation == "weighted_mean") {
    intensities <- intensities - peptide_level
  }

  # Abundance matrix: rows = groups, columns = samples
  n_groups <- ncol(membership)
  n_samples <- ncol(intensities)
  abundance <- matrix(
    NA_real_,
    nrow = n_groups,
    ncol = n_samples,
    dimnames = list(colnames(membership), colnames(intensities))
  )

  for (group in seq_len(n_groups)) {
    weights <- membership[, group]

    if (aggregation == "weighted_sum") {
      abundance[group, ] <- colSums(intensities * weights, na.rm = TRUE)
    } else {
      # Weighted mean over the peptides measured in each sample
      weighted_intensity_sum <- colSums(intensities * weights, na.rm = TRUE)
      weight_sum <- colSums((!is.na(intensities)) * weights, na.rm = TRUE)
      group_abundance <- weighted_intensity_sum / weight_sum
      group_abundance[weight_sum == 0] <- NA_real_
      if (center_peptides) {
        # Add back the group's level: the weighted mean of its peptides' levels
        group_abundance <- group_abundance +
          sum(weights * peptide_level, na.rm = TRUE) /
            sum(weights[is.finite(peptide_level)])
      }
      abundance[group, ] <- group_abundance
    }
  }

  group_labels <- rownames(abundance)
  if (is.null(group_labels)) {
    group_labels <- paste0("Group", seq_len(nrow(abundance)))
  }

  # Optionally average the replicate columns of each condition
  if (!is.null(condition_regex)) {
    column_conditions <- sub(condition_regex, "\\1", colnames(abundance), perl = TRUE)
    conditions <- unique(column_conditions)

    condition_abundance <- matrix(
      NA_real_,
      nrow = nrow(abundance),
      ncol = length(conditions),
      dimnames = list(group_labels, conditions)
    )
    for (condition in conditions) {
      condition_columns <- which(column_conditions == condition)
      condition_abundance[, condition] <- rowMeans(
        abundance[, condition_columns, drop = FALSE],
        na.rm = TRUE
      )
    }
    condition_abundance[is.nan(condition_abundance)] <- NA_real_

    abundance_wide <- data.table(
      ComplexoformGroup = group_labels,
      as.data.table(condition_abundance)
    )
  } else {
    abundance_wide <- data.table(
      ComplexoformGroup = group_labels,
      as.data.table(abundance)
    )
  }

  # Per group: contributing peptides and their total membership
  group_stats <- data.table(
    ComplexoformGroup = colnames(membership),
    NPeptides = colSums(membership > 0, na.rm = TRUE),
    MembershipWeightSum = colSums(membership, na.rm = TRUE)
  )

  list(
    abundance_wide = abundance_wide,
    group_stats = group_stats
  )
}


#' Quantify complexoform abundance from a VSClust result object
#'
#' Calls `quantify_complexoform_abundance()` with
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
  intensity_columns = NULL,
  id_columns = c("Gene name", "Peptide"),
  aggregation = c("weighted_mean", "weighted_sum"),
  membership_threshold = 0.5,
  normalize_membership = FALSE,
  center_peptides = TRUE,
  condition_regex = NULL
) {
  membership_matrix <- vsclust_result$ClustOut$Bestcl$membership
  peptide_data <- vsclust_result$original_data

  if (is.null(membership_matrix)) {
    stop("vsclust_result$ClustOut$Bestcl$membership is missing.")
  }
  if (is.null(peptide_data)) {
    stop("vsclust_result$original_data is missing.")
  }

  quantify_complexoform_abundance(
    peptide_data = peptide_data,
    membership_matrix = membership_matrix,
    intensity_columns = intensity_columns,
    id_columns = id_columns,
    aggregation = aggregation,
    membership_threshold = membership_threshold,
    normalize_membership = normalize_membership,
    center_peptides = center_peptides,
    condition_regex = condition_regex
  )
}
