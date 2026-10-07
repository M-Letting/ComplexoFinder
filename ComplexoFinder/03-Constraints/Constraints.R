library(data.table)


#' Core helper for peptide cannot-link constraints
#'
#' Implements the shared parsing, overlap filtering, and condition-mean logic
#' used by the various cannot-link constructors. The caller supplies the
#' scoring function and the rule that decides when to mark a pair as
#' cannot-link.
#'
#' @param data A `data.table` (or coercible) containing identifiers and
#'   intensity columns.
#' @param complex_name Character scalar. Name of the complex used to resolve
#'   multi-accession peptides (when `Position` contains `";"`).
#' @param id_cols Character vector. Columns in `data` that are pasted together
#'   (with `"_"`) to form `Identifier`.
#' @param identifier_info Character. Either a character vector of field names
#'   or a single underscore-delimited string. These names define the columns
#'   created when splitting `Identifier` by `"_"`.
#' @param intensity_cols Character vector of intensity columns. If `NULL`, uses
#'   all columns not in `id_cols`.
#' @param n_rep Integer. Number of replicates per condition.
#' @param n_cond Integer. Number of conditions. If `NULL`, inferred as
#'   `length(intensity_cols) / n_rep`.
#' @param min_shared_cond Integer. Minimum number of conditions with data in
#'   both peptides. Defaults to `floor(n_cond / 3)`.
#' @param cannot_link_th Numeric. Threshold passed to `cannot_link_rule`.
#' @param complex_info `data.table` with at least `complex_name` and
#'   `uniprot_id` columns, used to disambiguate multi-accession peptides.
#' @param verbose Logical. Whether to emit progress messages.
#' @param score_fun Function that maps two numeric vectors (condition means) to
#'   a scalar score.
#' @param cannot_link_rule Function with signature `(score, th)` returning
#'   `TRUE` when a cannot-link should be applied.
#' @param score_label Character. Label used in verbose messages.
#'
#' @return A logical matrix (`TRUE` = cannot-link) with dimnames set to the
#'   peptide identifiers.
#'
#' @keywords internal
create_cannot_link_core <- function(
  data,
  complex_name = NULL,
  id_cols = c("Gene name", "Peptide"),
  identifier_info = c(
    "Gene",
    "Accession",
    "Position",
    "Modification",
    "Datatype"
  ),
  intensity_cols = NULL,
  n_rep = 3,
  n_cond = NULL,
  min_shared_cond = NULL,
  cannot_link_th,
  complex_info = EBI_LT,
  verbose = FALSE,
  score_fun,
  cannot_link_rule,
  score_label = "score"
) {
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  missing_id_cols <- setdiff(id_cols, names(data))
  if (length(missing_id_cols) > 0) {
    stop(
      "Missing identifier columns: ",
      paste(missing_id_cols, collapse = ", ")
    )
  }
  if (is.null(intensity_cols)) {
    intensity_cols <- setdiff(names(data), id_cols)
  }
  if (length(intensity_cols) == 0) {
    stop("No intensity columns found after excluding id_cols.")
  }
  if (is.null(n_cond)) {
    n_cond <- length(intensity_cols) / n_rep
    if (n_cond != as.integer(n_cond)) {
      stop(
        "Cannot infer n_cond: length(intensity_cols) is not divisible by n_rep"
      )
    }
    n_cond <- as.integer(n_cond)
  }
  if (n_rep * n_cond != length(intensity_cols)) {
    stop("Number of intensity columns does not match n_rep * n_cond.")
  }
  if (is.null(min_shared_cond)) {
    min_shared_cond <- floor(n_cond / 3)
  }

  pick_col <- function(candidates) {
    direct_hit <- match(candidates, names(data))
    if (any(!is.na(direct_hit))) {
      return(candidates[which(!is.na(direct_hit))[1]])
    }
    lower_names <- tolower(names(data))
    lower_candidates <- tolower(candidates)
    lower_hit <- match(lower_candidates, lower_names)
    if (any(!is.na(lower_hit))) {
      return(names(data)[lower_hit[which(!is.na(lower_hit))[1]]])
    }
    NA_character_
  }

  grouping_col_data <- pick_col(c("Accession", "accession", "acc"))
  if (is.na(grouping_col_data)) {
    grouping_col_data <- pick_col(c("Gene name", "gene name"))
  }
  if (is.na(grouping_col_data)) {
    grouping_col_data <- pick_col(c("Gene", "gene"))
  }

  data[, Identifier := do.call(paste, c(.SD, sep = "_")), .SDcols = id_cols]
  peptide_ids <- sort(data$Identifier)

  identifier_cols <- identifier_info
  if (
    length(identifier_cols) == 1 && grepl("_", identifier_cols, fixed = TRUE)
  ) {
    identifier_cols <- strsplit(identifier_cols, "_", fixed = TRUE)[[1]]
  }
  if (length(identifier_cols) == 0) {
    stop("identifier_info must define at least one field name")
  }

  split_vals <- tstrsplit(peptide_ids, "_", fixed = TRUE, fill = NA_character_)
  if (length(split_vals) < length(identifier_cols)) {
    split_vals <- c(
      split_vals,
      rep(
        list(rep(NA_character_, length(peptide_ids))),
        length(identifier_cols) - length(split_vals)
      )
    )
  }
  split_vals <- split_vals[seq_len(length(identifier_cols))]

  peptide_info <- data.table(Identifier = peptide_ids)
  for (k in seq_along(identifier_cols)) {
    peptide_info[, (identifier_cols[k]) := split_vals[[k]]]
  }

  if (!is.na(grouping_col_data)) {
    peptide_info[,
      GroupKey := data[match(peptide_ids, Identifier), get(grouping_col_data)]
    ]
  } else {
    peptide_info[, GroupKey := NA_character_]
  }

  peptide_info[, `:=`(
    Start = NA_integer_,
    End = NA_integer_
  )]

  grouping_label <- if (!is.na(grouping_col_data)) {
    grouping_col_data
  } else {
    "<none>"
  }
  accession_col <- if ("Accession" %in% names(peptide_info)) {
    "Accession"
  } else if (length(identifier_cols) >= 2) {
    identifier_cols[2]
  } else {
    NA_character_
  }
  position_col <- if ("Position" %in% names(peptide_info)) {
    "Position"
  } else {
    pos_hits <- grep(
      "^position$|pos",
      names(peptide_info),
      ignore.case = TRUE,
      value = TRUE
    )
    if (length(pos_hits) > 0) pos_hits[1] else NA_character_
  }

  if (verbose) {
    message("Parsing peptide positions...")
  }

  for (i in seq_len(nrow(peptide_info))) {
    if (is.na(position_col)) {
      break
    }
    position_str <- peptide_info[[position_col]][i]
    position_str <- as.character(position_str)[1]
    if (!isTRUE(nzchar(position_str))) {
      next
    }

    if (grepl(";|\\|", position_str)) {
      if (verbose) {
        message("  Multi-accession peptide: ", position_str)
      }
      if (is.na(accession_col)) {
        next
      }
      accession_str <- peptide_info[[accession_col]][i]
      accession_str <- as.character(accession_str)[1]
      if (!isTRUE(nzchar(accession_str))) {
        next
      }

      accessions <- unlist(strsplit(accession_str, ";|\\|"))
      accessions <- trimws(accessions)

      complex_info_subset <- complex_info[
        complex_info$complex_name == complex_name,
      ]
      # Normalize isoform-1 to canonical before matching (Q13936-1 == Q13936)
      norm_acc <- function(a) sub("-1$", "", a)
      in_complex <- norm_acc(accessions) %in%
        norm_acc(complex_info_subset$uniprot_id)

      if (verbose) {
        message(
          "    Accessions in complex: ",
          paste(accessions[in_complex], collapse = ", ")
        )
      }

      if (sum(in_complex) == 1) {
        complex_idx <- which(in_complex)[1]
        position_entries <- unlist(strsplit(position_str, ";|\\|"))
        position_entries <- trimws(position_entries)
        position_entry <- position_entries[complex_idx]

        m <- regmatches(
          position_entry,
          regexec("\\[(\\d+)\\s*[-–]\\s*(\\d+)\\]", position_entry, perl = TRUE)
        )
        if (length(m[[1]]) < 3) {
          m <- regmatches(
            position_entry,
            regexec("(\\d+)\\s*[-–]\\s*(\\d+)", position_entry, perl = TRUE)
          )
        }

        if (length(m) > 0 && length(m[[1]]) >= 3) {
          peptide_info$Start[i] <- as.integer(m[[1]][2])
          peptide_info$End[i] <- as.integer(m[[1]][3])
        }
      }
    } else {
      m <- regmatches(
        position_str,
        regexec("\\[(\\d+)\\s*[-–]\\s*(\\d+)\\]", position_str, perl = TRUE)
      )
      if (length(m[[1]]) < 3) {
        m <- regmatches(
          position_str,
          regexec("(\\d+)\\s*[-–]\\s*(\\d+)", position_str, perl = TRUE)
        )
      }
      if (length(m) > 0 && length(m[[1]]) >= 3) {
        peptide_info$Start[i] <- as.integer(m[[1]][2])
        peptide_info$End[i] <- as.integer(m[[1]][3])
      }
    }
  }

  if (verbose) {
    message(
      "Parsed ",
      sum(!is.na(peptide_info$Start)),
      " of ",
      nrow(peptide_info),
      " peptide positions"
    )
  }

  n_peptides <- length(peptide_ids)
  constraints <- matrix(FALSE, nrow = n_peptides, ncol = n_peptides)
  rownames(constraints) <- peptide_ids
  colnames(constraints) <- peptide_ids

  data_matrix <- as.matrix(data[
    match(peptide_ids, data$Identifier),
    ..intensity_cols
  ])
  rownames(data_matrix) <- peptide_ids

  valid_positions <- !is.na(peptide_info$Start) & !is.na(peptide_info$End)

  if (verbose) {
    message(
      "Building constraints for ",
      sum(valid_positions),
      " peptides with valid positions..."
    )
  }

  if (verbose && !is.na(grouping_col_data)) {
    message("Grouping comparisons by: ", grouping_label)
  }

  group_vec <- peptide_info$GroupKey
  if (is.na(grouping_col_data)) {
    group_vec <- rep("ALL", nrow(peptide_info))
  }

  unique_groups <- unique(group_vec)

  for (group in unique_groups) {
    if (is.na(group)) {
      next
    }
    group_indices <- which(group_vec == group & valid_positions)
    if (length(group_indices) < 2) {
      next
    }

    if (verbose) {
      message(
        "  Processing ",
        length(group_indices),
        " peptides for group: ",
        group
      )
    }

    for (idx_i in seq_along(group_indices)) {
      for (idx_j in seq_along(group_indices)) {
        i <- group_indices[idx_i]
        j <- group_indices[idx_j]

        if (i >= j) {
          next
        }

        start_i <- peptide_info$Start[i]
        end_i <- peptide_info$End[i]
        start_j <- peptide_info$Start[j]
        end_j <- peptide_info$End[j]

        no_overlap <- (end_i < start_j) || (end_j < start_i)
        if (no_overlap) {
          next
        }

        intensities_i <- data_matrix[i, ]
        intensities_j <- data_matrix[j, ]

        cond_matrix_i <- matrix(
          intensities_i,
          nrow = n_rep,
          ncol = n_cond,
          byrow = FALSE
        )
        cond_matrix_j <- matrix(
          intensities_j,
          nrow = n_rep,
          ncol = n_cond,
          byrow = FALSE
        )

        cond_mean_i <- colMeans(cond_matrix_i, na.rm = TRUE)
        cond_mean_j <- colMeans(cond_matrix_j, na.rm = TRUE)

        shared_conditions <- sum(!is.na(cond_mean_i) & !is.na(cond_mean_j))
        if (shared_conditions < min_shared_cond) {
          next
        }

        overlap_start <- max(start_i, start_j)
        overlap_end <- min(end_i, end_j)
        overlap_length <- overlap_end - overlap_start + 1

        len_i <- end_i - start_i + 1
        len_j <- end_j - start_j + 1
        overlap_fraction <- overlap_length / min(len_i, len_j)
        if (overlap_fraction <= 0.5) {
          next
        }

        score <- score_fun(cond_mean_i, cond_mean_j)

        if (!is.na(score) && cannot_link_rule(score, cannot_link_th)) {
          constraints[i, j] <- TRUE
          constraints[j, i] <- TRUE
        }
      }
    }
  }

  if (verbose) {
    message(
      "Created ",
      sum(constraints) / 2,
      " cannot-link constraints using ",
      score_label,
      " threshold ",
      cannot_link_th
    )
  }

  constraints
}


#' Create peptide cannot-link constraints from overlap and correlation
#'
#' Constructs a peptide × peptide cannot-link matrix based on (i) overlapping
#' peptide positions within the same gene/protein and (ii) low correlation of
#' condition means.
#'
#' The function first creates a single string identifier by pasting the columns
#' listed in `id_cols` with `"_"`. It then splits that identifier into the
#' fields specified by `identifier_info` and adds `Start`/`End` columns parsed
#' from the `Position` field.
#'
#' @param data A `data.table` (or coercible) containing identifiers and
#'   intensity columns.
#' @param complex_name Character scalar. Name of the complex used to resolve
#'   multi-accession peptides (when `Position` contains `";"`).
#' @param id_cols Character vector. Columns in `data` that are pasted together
#'   (with `"_"`) to form `Identifier`.
#' @param identifier_info Character. Either a character vector of field names
#'   (e.g. `c("Gene","Accession","Position","Modification","Datatype")`) or a
#'   single underscore-delimited string (e.g.
#'   `"Gene_Accession_Position_Modification_Datatype"`). These names define the
#'   columns created when splitting `Identifier` by `"_"`.
#' @param intensity_cols Character vector of intensity columns. If `NULL`, uses
#'   all columns not in `id_cols`.
#' @param n_rep Integer. Number of replicates per condition.
#' @param n_cond Integer. Number of conditions. If `NULL`, inferred as
#'   `length(intensity_cols) / n_rep`.
#' @param min_shared_cond Integer. Minimum number of conditions with data in
#'   both peptides. Defaults to `floor(n_cond / 3)`.
#' @param cannot_link_th Numeric. Apply cannot-link when correlation is below
#'   this threshold.
#' @param complex_info `data.table` with at least `complex_name` and
#'   `uniprot_id` columns, used to disambiguate multi-accession peptides.
#' @param verbose Logical. Whether to emit progress messages.
#'
#' @return A logical matrix (`TRUE` = cannot-link) with dimnames set to the
#'   peptide identifiers.
#'
#' @export
create_cannot_link_Pearson <- function(
  data,
  complex_name = NULL,
  id_cols = c("Gene name", "Peptide"),
  identifier_info = c(
    "Gene",
    "Accession",
    "Position",
    "Modification",
    "Datatype"
  ),
  intensity_cols = NULL,
  n_rep = 3,
  n_cond = NULL,
  min_shared_cond = NULL,
  cannot_link_th = 0.6,
  complex_info = EBI_LT,
  verbose = FALSE
) {
  score_fun <- function(x, y) {
    valid <- is.finite(x) & is.finite(y)
    if (sum(valid) < 3) {
      return(NA_real_)
    }
    cor(x[valid], y[valid])
  }

  create_cannot_link_core(
    data = data,
    complex_name = complex_name,
    id_cols = id_cols,
    identifier_info = identifier_info,
    intensity_cols = intensity_cols,
    n_rep = n_rep,
    n_cond = n_cond,
    min_shared_cond = min_shared_cond,
    cannot_link_th = cannot_link_th,
    complex_info = complex_info,
    verbose = verbose,
    score_fun = score_fun,
    cannot_link_rule = function(score, th) score < th,
    score_label = "Pearson"
  )
}


#' Create peptide cannot-link constraints using robust z-scored RMSD
#'
#' Constructs a peptide × peptide cannot-link matrix using the same overlap and
#' location logic as `create_cannot_link()`, but replaces Pearson correlation
#' with robust shape discordance measured as RMSD between robust-z transformed
#' condition means.
#'
#' @inheritParams create_cannot_link_Pearson
#' @param cannot_link_th Numeric. Apply cannot-link when robust-z RMSD is above
#'   this threshold.
#'
#' @return A logical matrix (`TRUE` = cannot-link) with dimnames set to the
#'   peptide identifiers.
#'
#' @export
create_cannot_link_ZRMSD <- function(
  data,
  complex_name = NULL,
  id_cols = c("Gene name", "Peptide"),
  identifier_info = c(
    "Gene",
    "Accession",
    "Position",
    "Modification",
    "Datatype"
  ),
  intensity_cols = NULL,
  n_rep = 3,
  n_cond = NULL,
  min_shared_cond = NULL,
  cannot_link_th = 1.0,
  complex_info = EBI_LT,
  verbose = FALSE
) {
  robust_z <- function(x) {
    med <- median(x, na.rm = TRUE)
    mad_val <- mad(x, constant = 1, na.rm = TRUE)
    if (!is.finite(mad_val) || mad_val == 0) {
      return(rep(NA_real_, length(x)))
    }
    (x - med) / mad_val
  }

  score_fun <- function(x, y) {
    valid <- is.finite(x) & is.finite(y)
    if (sum(valid) < 3) {
      return(NA_real_)
    }

    z_x <- robust_z(x)
    z_y <- robust_z(y)

    valid <- valid & is.finite(z_x) & is.finite(z_y)
    if (sum(valid) < 3) {
      return(NA_real_)
    }

    sqrt(mean((z_x[valid] - z_y[valid])^2))
  }

  create_cannot_link_core(
    data = data,
    complex_name = complex_name,
    id_cols = id_cols,
    identifier_info = identifier_info,
    intensity_cols = intensity_cols,
    n_rep = n_rep,
    n_cond = n_cond,
    min_shared_cond = min_shared_cond,
    cannot_link_th = cannot_link_th,
    complex_info = complex_info,
    verbose = verbose,
    score_fun = score_fun,
    cannot_link_rule = function(score, th) score > th,
    score_label = "robust z-RMSD"
  )
}


#' Create peptide cannot-link constraints using Concordance Correlation Coefficient
#'
#' Constructs a peptide × peptide cannot-link matrix using the same overlap and
#' location logic as `create_cannot_link()`, but replaces Pearson correlation
#' with Lin's concordance correlation coefficient (CCC), which penalizes both
#' bias and scale differences.
#'
#' @inheritParams create_cannot_link_Pearson
#' @param cannot_link_th Numeric. Apply cannot-link when CCC is below this
#'   threshold.
#'
#' @return A logical matrix (`TRUE` = cannot-link) with dimnames set to the
#'   peptide identifiers.
#'
#' @export
create_cannot_link_CCC <- function(
  data,
  complex_name = NULL,
  id_cols = c("Gene name", "Peptide"),
  identifier_info = c(
    "Gene",
    "Accession",
    "Position",
    "Modification",
    "Datatype"
  ),
  intensity_cols = NULL,
  n_rep = 3,
  n_cond = NULL,
  min_shared_cond = NULL,
  cannot_link_th = 0.5,
  complex_info = EBI_LT,
  verbose = FALSE
) {
  score_fun <- function(x, y) {
    valid <- is.finite(x) & is.finite(y)
    if (sum(valid) < 3) {
      return(NA_real_)
    }

    x <- x[valid]
    y <- y[valid]

    mean_x <- mean(x)
    mean_y <- mean(y)
    var_x <- var(x)
    var_y <- var(y)
    cov_xy <- cov(x, y)

    denom <- var_x + var_y + (mean_x - mean_y)^2
    if (!is.finite(denom) || denom == 0) {
      return(NA_real_)
    }

    (2 * cov_xy) / denom
  }

  create_cannot_link_core(
    data = data,
    complex_name = complex_name,
    id_cols = id_cols,
    identifier_info = identifier_info,
    intensity_cols = intensity_cols,
    n_rep = n_rep,
    n_cond = n_cond,
    min_shared_cond = min_shared_cond,
    cannot_link_th = cannot_link_th,
    complex_info = complex_info,
    verbose = verbose,
    score_fun = score_fun,
    cannot_link_rule = function(score, th) score < th,
    score_label = "CCC"
  )
}


#' Create VSClust Restriction Matrix from dCF Assignments
#'
#' Converts dCF cluster assignments and cannot-link constraints into a
#' peptide × cluster restriction matrix for semi-supervised VSClust clustering.
#'
#' @param data data.table with dCF assignments and peptide identifiers
#' @param cannotlink_matrix Peptide × peptide cannot-link constraint matrix (from create_peptide_constraints)
#' @param dCF_col Character. Column name containing dCF assignments (default: "dCF")
#' @param id_cols Character vector. Columns to create peptide identifiers (default: c("Gene name", "Peptide"))
#' @param verbose Logical. Whether to print progress messages
#'
#' @return Matrix with peptides as rows, clusters as columns
#'   - Values are TRUE or FALSE
#'   - FALSE = peptide is allowed in this cluster
#'   - TRUE = peptide is forbidden from this cluster
#'   - Rownames = peptide identifiers (Gene_Peptide format)
#'   - Colnames = cluster numbers (1, 2, 3, ..., n_clusters)
#'   - dCF0 (canonical) is included as cluster 1
#'   - dCF1, dCF2, ... become clusters 2, 3, ...
#'
#' @details
#' Algorithm:
#' 1. Initialize restriction matrix with all FALSE (all peptides allowed in all clusters)
#' 2. For each cannot-link pair (i,j):
#'    - Set restriction[i, cluster_j] = TRUE (i forbidden from j's cluster)
#'    - Set restriction[j, cluster_i] = TRUE (j forbidden from i's cluster)
#' 3. Handle conflicts: if peptide forbidden from all clusters, allow into all
#'
#' @export
dCF_to_constraints <- function(
  data,
  cannotlink_matrix,
  dCF_col = "dCF",
  id_cols = c("Gene name", "Peptide"),
  verbose = FALSE
) {
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }

  # Check required columns
  if (!dCF_col %in% names(data)) {
    stop("Column '", dCF_col, "' not found in data")
  }

  missing_cols <- setdiff(id_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing identifier columns: ", paste(missing_cols, collapse = ", "))
  }

  # Create unique peptide identifier
  data[, Identifier := do.call(paste, c(.SD, sep = "_")), .SDcols = id_cols]

  # Filter out singletons (dCF-1) ONLY - keep dCF0 as canonical cluster
  data_filtered <- data[get(dCF_col) != "dCF-1"]

  if (nrow(data_filtered) == 0) {
    stop("No valid dCF assignments found (all peptides are dCF-1)")
  }

  # Extract dCF numbers (dCF0 = 0, dCF1 = 1, etc.)
  data_filtered[, dCF_num := as.integer(gsub("dCF", "", get(dCF_col)))]

  # Get unique clusters INCLUDING dCF0 (0, 1, 2, 3, ...)
  unique_clusters <- sort(unique(data_filtered$dCF_num))
  n_clusters <- length(unique_clusters)

  if (n_clusters == 0) {
    stop("No valid clusters found (all peptides are dCF-1)")
  }

  # Get peptides from data (match cannotlink_matrix rownames)
  peptides_in_data <- data_filtered$Identifier
  peptides_in_constraints <- rownames(cannotlink_matrix)
  all_peptides <- intersect(peptides_in_constraints, peptides_in_data)
  n_peptides <- length(all_peptides)

  if (verbose) {
    message("Creating VSClust restriction matrix:")
    message("  ", n_peptides, " peptides")
    message(
      "  ",
      n_clusters,
      " clusters (dCF",
      paste(unique_clusters, collapse = ", dCF"),
      ")"
    )
  }

  # Create lookup for peptide -> cluster assignment
  pseudo_assignment <- setNames(
    data_filtered[match(all_peptides, Identifier), dCF_num],
    all_peptides
  )

  # Initialize restriction matrix: all FALSE (all peptides allowed everywhere)
  restriction_matrix <- matrix(
    FALSE,
    nrow = n_peptides,
    ncol = n_clusters
  )
  rownames(restriction_matrix) <- all_peptides
  colnames(restriction_matrix) <- as.character(seq_len(n_clusters))

  if (verbose) {
    for (i in seq_along(unique_clusters)) {
      cluster_num <- unique_clusters[i]
      n_in_cluster <- sum(pseudo_assignment == cluster_num, na.rm = TRUE)
      message(
        "  dCF",
        cluster_num,
        " -> Cluster ",
        i,
        ": ",
        n_in_cluster,
        " peptides"
      )
    }
  }

  # Iterate over all peptide pairs and enforce cannot-link constraints
  if (!is.null(cannotlink_matrix) && sum(cannotlink_matrix) > 0) {
    if (verbose) {
      message(
        "  Enforcing ",
        sum(cannotlink_matrix) / 2,
        " cannot-link constraints..."
      )
    }

    n_restrictions_added <- 0

    # Only check peptides that are in both matrices
    for (i in seq_along(all_peptides)) {
      peptide_i <- all_peptides[i]
      cluster_i <- pseudo_assignment[peptide_i]

      # Skip if peptide i is not assigned to a cluster (NA)
      if (is.na(cluster_i)) {
        next
      }

      # Find which cluster index this corresponds to (dCF0 -> index 1, dCF1 -> index 2, etc.)
      cluster_i_idx <- which(unique_clusters == cluster_i)
      if (length(cluster_i_idx) == 0) {
        next
      }

      for (j in seq(i + 1, length(all_peptides))) {
        if (j > length(all_peptides)) {
          break
        }

        peptide_j <- all_peptides[j]

        # Check if there's a cannot-link constraint
        if (cannotlink_matrix[peptide_i, peptide_j]) {
          cluster_j <- pseudo_assignment[peptide_j]

          # Skip if peptide j is not assigned to a cluster (NA)
          if (is.na(cluster_j)) {
            next
          }

          # Find which cluster index this corresponds to
          cluster_j_idx <- which(unique_clusters == cluster_j)
          if (length(cluster_j_idx) == 0) {
            next
          }

          # Set restrictions
          restriction_matrix[peptide_i, cluster_j_idx] <- TRUE
          restriction_matrix[peptide_j, cluster_i_idx] <- TRUE
          n_restrictions_added <- n_restrictions_added + 2
        }
      }
    }

    if (verbose) {
      message("  Added ", n_restrictions_added, " restriction entries")
    }
  }

  # Handle conflicts: peptides forbidden from all clusters
  all_forbidden <- rowSums(restriction_matrix) == n_clusters
  if (any(all_forbidden)) {
    if (verbose) {
      message(
        "  WARNING: ",
        sum(all_forbidden),
        " peptides forbidden from all clusters - allowing into all"
      )
    }
    restriction_matrix[all_forbidden, ] <- FALSE
  }

  # Report statistics
  if (verbose) {
    cluster_sizes <- colSums(!restriction_matrix)
    message(
      "  Allowed peptides per cluster: ",
      paste(cluster_sizes, collapse = ", ")
    )
  }

  restriction_matrix
}


#' Create VSClust Restriction Matrix from Initial Clustering Results
#'
#' Runs an initial unconstrained VSClust clustering, then uses the clustering
#' assignments combined with cannot-link constraints to build a restriction matrix.
#' This restriction matrix can then be used for a second, constrained clustering pass.
#'
#' @param data data.table with intensity values and peptide identifiers
#' @param cannotlink_matrix Peptide × peptide cannot-link constraint matrix (from create_peptide_constraints)
#' @param n_clusters Integer. Number of clusters to find (typically number of dCF excluding dCF-1)
#' @param n_rep Integer. Number of replicates per condition
#' @param n_cond Integer. Number of conditions
#' @param id_cols Character vector. Columns to create peptide identifiers (default: c("Gene name", "Peptide"))
#' @param value_cols Character vector. Column names containing intensity values
#' @param verbose Logical. Whether to print progress messages
#'
#' @return List with two elements:
#'   - restriction_matrix: Matrix with peptides as rows, clusters as columns
#'     * Values are TRUE (forbidden) or FALSE (allowed)
#'     * Rownames = peptide identifiers (Gene_Peptide format)
#'     * Colnames = cluster numbers (1, 2, 3, ..., n_clusters)
#'   - initial_clustering: The unconstrained VSClust clustering results
#'
#' @details
#' Algorithm:
#' 1. Run unconstrained VSClust clustering with n_clusters
#' 2. Initialize restriction matrix with all FALSE (all peptides allowed everywhere)
#' 3. For each cannot-link pair (pep1, pep2):
#'    a. If assigned to different clusters (cluster1, cluster2):
#'       - Set restriction[pep1, cluster2] = TRUE (pep1 forbidden from cluster2)
#'       - Set restriction[pep2, cluster1] = TRUE (pep2 forbidden from cluster1)
#'    b. If assigned to same cluster (conflict - unlikely but possible):
#'       - Compare membership values for both peptides in that cluster
#'       - Keep peptide with higher membership, forbid the other
#'       - Set restriction[lower_membership_peptide, shared_cluster] = TRUE
#'
#' @export
vsclust_to_restrictions <- function(
  data,
  cannotlink_matrix,
  n_clusters,
  n_rep = 3,
  n_cond = NULL,
  id_cols = c("Gene name", "Peptide"),
  value_cols = NULL,
  verbose = FALSE
) {
  if (!requireNamespace("vsclust", quietly = TRUE)) {
    stop(
      "Package 'vsclust' is required. Install with: BiocManager::install('vsclust')"
    )
  }

  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }

  # Check required columns
  missing_id_cols <- setdiff(id_cols, names(data))
  if (length(missing_id_cols) > 0) {
    stop(
      "Missing identifier columns: ",
      paste(missing_id_cols, collapse = ", ")
    )
  }

  if (is.null(value_cols)) {
    # Auto-detect intensity columns (exclude id_cols)
    value_cols <- setdiff(names(data), id_cols)
  }

  if (is.null(n_cond)) {
    n_cond <- length(value_cols) / n_rep
  }

  if (verbose) {
    message("Step 1: Running initial unconstrained VSClust clustering...")
    message("  Peptides: ", nrow(data))
    message("  Clusters: ", n_clusters)
    message("  Conditions: ", n_cond, " x ", n_rep, " replicates")
  }

  # Run initial unconstrained clustering
  initial_clustering <- vsclust_on_complex(
    data = data,
    id_cols = id_cols,
    intensity_cols = value_cols,
    n_rep = n_rep,
    n_cond = n_cond,
    n_clusters = n_clusters,
    grouped_replicates = TRUE,
    restriction_matrix = NULL # No constraints for initial run
  )

  # Extract cluster assignments and membership values
  cluster_assignments <- initial_clustering$ClustOut$Bestcl$cluster
  membership_matrix <- initial_clustering$ClustOut$Bestcl$membership

  # Create peptide identifiers to match cannotlink_matrix rownames
  data[, peptide_id := paste(get(id_cols[1]), get(id_cols[2]), sep = "_")]
  peptide_ids <- data$peptide_id
  names(cluster_assignments) <- peptide_ids
  rownames(membership_matrix) <- peptide_ids

  if (verbose) {
    message(
      "Step 2: Building restriction matrix from clustering + cannot-link..."
    )
    message(
      "  Cluster sizes: ",
      paste(table(cluster_assignments), collapse = ", ")
    )
  }

  # Initialize restriction matrix: all FALSE (all allowed)
  n_peptides <- length(peptide_ids)
  restriction_matrix <- matrix(
    FALSE,
    nrow = n_peptides,
    ncol = n_clusters,
    dimnames = list(peptide_ids, 1:n_clusters)
  )

  # Apply cannot-link constraints
  n_restrictions_added <- 0
  n_conflicts_resolved <- 0

  # Find all cannot-link pairs
  cannotlink_pairs <- which(cannotlink_matrix, arr.ind = TRUE)
  # Only keep upper triangle to avoid processing each pair twice
  cannotlink_pairs <- cannotlink_pairs[
    cannotlink_pairs[, 1] < cannotlink_pairs[, 2],
    ,
    drop = FALSE
  ]

  if (nrow(cannotlink_pairs) == 0) {
    if (verbose) {
      message("  No cannot-link constraints to apply")
    }
  } else {
    if (verbose) {
      message("  Processing ", nrow(cannotlink_pairs), " cannot-link pairs...")
    }

    for (i in seq_len(nrow(cannotlink_pairs))) {
      row_idx <- cannotlink_pairs[i, 1]
      col_idx <- cannotlink_pairs[i, 2]

      peptide_i <- rownames(cannotlink_matrix)[row_idx]
      peptide_j <- rownames(cannotlink_matrix)[col_idx]

      # Skip if peptides not in clustering results (filtered out)
      if (
        !peptide_i %in% names(cluster_assignments) ||
          !peptide_j %in% names(cluster_assignments)
      ) {
        next
      }

      cluster_i <- cluster_assignments[peptide_i]
      cluster_j <- cluster_assignments[peptide_j]

      if (cluster_i != cluster_j) {
        # Case A: Peptides in different clusters
        # Forbid each peptide from the other's cluster
        restriction_matrix[peptide_i, cluster_j] <- TRUE
        restriction_matrix[peptide_j, cluster_i] <- TRUE
        n_restrictions_added <- n_restrictions_added + 2
      } else {
        # Case B: Peptides in same cluster (conflict!)
        # Allow only the peptide with higher membership value
        membership_i <- membership_matrix[peptide_i, cluster_i]
        membership_j <- membership_matrix[peptide_j, cluster_i]

        if (membership_i > membership_j) {
          # Forbid peptide_j from this cluster
          restriction_matrix[peptide_j, cluster_i] <- TRUE
          if (verbose) {
            message(
              "  Conflict: ",
              peptide_i,
              " (mem=",
              round(membership_i, 3),
              ") vs ",
              peptide_j,
              " (mem=",
              round(membership_j, 3),
              ") in cluster ",
              cluster_i,
              " -> forbid ",
              peptide_j
            )
          }
        } else {
          # Forbid peptide_i from this cluster
          restriction_matrix[peptide_i, cluster_i] <- TRUE
          if (verbose) {
            message(
              "  Conflict: ",
              peptide_j,
              " (mem=",
              round(membership_j, 3),
              ") vs ",
              peptide_i,
              " (mem=",
              round(membership_i, 3),
              ") in cluster ",
              cluster_i,
              " -> forbid ",
              peptide_i
            )
          }
        }
        n_restrictions_added <- n_restrictions_added + 1
        n_conflicts_resolved <- n_conflicts_resolved + 1
      }
    }
  }

  # Handle edge case: peptides forbidden from all clusters
  all_forbidden <- rowSums(restriction_matrix) == n_clusters
  if (any(all_forbidden)) {
    if (verbose) {
      message(
        "  WARNING: ",
        sum(all_forbidden),
        " peptides forbidden from all clusters - allowing into all"
      )
    }
    restriction_matrix[all_forbidden, ] <- FALSE
  }

  # Report statistics
  if (verbose) {
    message("  Added ", n_restrictions_added, " restriction entries")
    if (n_conflicts_resolved > 0) {
      message("  Resolved ", n_conflicts_resolved, " same-cluster conflicts")
    }
    cluster_sizes <- colSums(!restriction_matrix)
    message(
      "  Allowed peptides per cluster: ",
      paste(cluster_sizes, collapse = ", ")
    )
  }

  return(list(
    restriction_matrix = restriction_matrix,
    initial_clustering = initial_clustering
  ))
}
