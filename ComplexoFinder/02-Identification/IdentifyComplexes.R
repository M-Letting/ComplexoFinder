library(data.table)

#' Find Protein Complexes in Dataset
#'
#' Identifies protein complexes in a list of datasets by matching protein IDs
#' against a lookup table (e.g., EBI Complex Portal database).
#'
#' Matching is accession-based with isoform awareness. UniProt isoform-1
#' accessions in the lookup (e.g. \code{Q13936-1}) are treated as equivalent
#' to the canonical base accession (\code{Q13936}), because isoform 1 is the
#' canonical form and is routinely reported without a suffix by proteomics
#' software. Isoforms 2 and higher (e.g. \code{Q13936-2}) require an exact
#' match: a data row must carry the same isoform-specific accession to be
#' counted as a complex member.
#'
#' @param data_list Named list of data.tables containing proteomics data.
#' @param lookup_table data.table with complex annotations containing columns:
#'   \code{complex_id}, \code{complex_name}, and the protein identifier column
#'   specified by \code{id_col_lookup}. UniProt accessions with isoform
#'   suffixes (e.g. \code{Q13936-1}, \code{Q13936-2}) are supported.
#' @param id_col_data Character. Column name in data containing protein
#'   accessions (e.g. \code{"Accession"}). Multi-protein rows with
#'   semicolon-separated values are expanded before matching.
#' @param id_col_lookup Character. Column name in \code{lookup_table} for
#'   accession matching (e.g. \code{"uniprot_id"}).
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{complexes}: Named list of data.tables, one per complex,
#'       containing all matching rows from the input datasets.
#'     \item \code{summary}: data.frame with per-complex statistics:
#'       proteins found, total annotated members, fraction detected, and
#'       total peptide rows.
#'   }
#'
#' @examples
#' \dontrun{
#' result <- find_complexes(
#'   data_list = H9_data,
#'   lookup_table = EBI_LT,
#'   id_col_data = "Accession",
#'   id_col_lookup = "uniprot_id"
#' )
#' }
#'
#' @export
find_complexes <- function(
  data_list,
  lookup_table,
  id_col_data,
  id_col_lookup
) {
  # Remove duplicate entries in lookup table based id_col_lookup
  lookup_table <- lookup_table[
    !duplicated(paste(lookup_table$complex_id, lookup_table[[id_col_lookup]])),
  ]

  # Convert lookup to data.table for fast operations
  setDT(lookup_table)

  # Normalized match key: isoform-1 accessions (X-1) also match canonical (X)
  lookup_table[, .match_id := sub("-1$", "", get(id_col_lookup))]

  # Initialize result list organized by complex
  all_complexes <- list()

  # Process each dataset
  for (i in seq_along(data_list)) {
    dataset_name <- names(data_list)[i]
    data <- copy(data_list[[i]])
    setDT(data)

    # Add row index for tracking
    data[, .row_idx := .I]

    # Split multi-ID rows into separate rows for matching
    # Each row with "ID1; ID2" becomes multiple rows with single IDs
    expanded <- data[,
      .(
        protein_id = trimws(unlist(strsplit(get(id_col_data), ";"))),
        .row_idx = .row_idx
      ),
      by = .row_idx
    ]
    expanded[, .match_id := sub("-1$", "", protein_id)]

    # Join with lookup table to find complex memberships
    # This is a vectorized operation - much faster than nested loops
    matched <- lookup_table[
      expanded,
      on = ".match_id",
      nomatch = NULL,
      allow.cartesian = TRUE
    ]
    matched[, .match_id := NULL]

    # For each complex that has matches
    if (nrow(matched) > 0) {
      for (complex_id in unique(matched$complex_id)) {
        complex_name <- unique(matched$complex_name[
          matched$complex_id == complex_id
        ])[1]

        # Get original row indices that match this complex
        matching_row_indices <- unique(matched$.row_idx[
          matched$complex_id == complex_id
        ])

        # Extract original (unsplit) rows from data
        matched_data <- copy(data[matching_row_indices, ])
        matched_data[, .row_idx := NULL]
        matched_data$complex_id <- complex_id
        matched_data$complex_name <- complex_name

        # Initialize complex entry if it doesn't exist
        if (is.null(all_complexes[[complex_name]])) {
          all_complexes[[complex_name]] <- list()
        }

        # Add this dataset's proteins to the complex
        all_complexes[[complex_name]][[dataset_name]] <- matched_data
      }
    }
  }

  lookup_table[, .match_id := NULL]

  # Merge all datasets for each complex
  all_complexes <- lapply(all_complexes, function(complex_data) {
    do.call(rbind, complex_data)
  })

  # Remove complexes with only one unique protein/accession in total
  # Need to count unique protein IDs (handling multiple IDs per row)
  all_complexes <- Filter(
    function(df) {
      # Extract all unique protein IDs
      all_ids <- unique(unlist(lapply(df[[id_col_data]], function(x) {
        if (is.na(x) || x == "") {
          return(character(0))
        }
        trimws(strsplit(x, ";")[[1]])
      })))
      length(all_ids) > 1
    },
    all_complexes
  )

  # Create summary data frame
  summary_df <- data.frame(
    complex_name = names(all_complexes),
    proteins_found = vapply(
      all_complexes,
      function(df) {
        # Count unique protein IDs (handling multiple per row)
        all_ids <- unique(unlist(lapply(df[[id_col_data]], function(x) {
          if (is.na(x) || x == "") {
            return(character(0))
          }
          trimws(strsplit(x, ";")[[1]])
        })))
        length(all_ids)
      },
      integer(1)
    ),
    proteins_total = vapply(
      names(all_complexes),
      function(nm) {
        length(unique(sub("-1$", "", lookup_table[[id_col_lookup]][
          lookup_table$complex_name == nm &
            !is.na(lookup_table[[id_col_lookup]])
        ])))
      },
      integer(1)
    ),
    peptides_total = vapply(all_complexes, nrow, integer(1))
  )
  summary_df$fraction_found <- summary_df$proteins_found /
    summary_df$proteins_total

  # Return list with complexes and summary
  list(
    complexes = all_complexes,
    summary = summary_df
  )
}

#' Prune Complexes Using Statistical Enrichment
#'
#' Filters complexes using a hypergeometric enrichment test to retain only
#' those complexes whose detected members are significantly over-represented
#' relative to a proteome-wide background.
#'
#' Isoform handling mirrors \code{\link{find_complexes}}: isoform-1 accessions
#' in the lookup (e.g. \code{Q13936-1}) are normalized to the canonical base
#' accession (\code{Q13936}) before counting annotated members and computing
#' overlap with detected proteins. Isoforms 2+ are kept as-is and must match
#' exactly in the data. The \code{universe} should therefore contain canonical
#' (base) accessions without isoform suffixes, such as those obtained from
#' the UniProt reviewed human proteome (\code{human_proteome$Accession}).
#'
#' @param all_complexes Named list of complex data.tables, as returned by the
#'   \code{complexes} element of \code{\link{find_complexes}}.
#' @param lookup_table Complex lookup table (same object passed to
#'   \code{\link{find_complexes}}).
#' @param id_col_data Character. Column name for protein accessions in the
#'   complex data.tables (e.g. \code{"Accession"}).
#' @param id_col_lookup Character. Column name for protein accessions in
#'   \code{lookup_table} (e.g. \code{"uniprot_id"}).
#' @param universe Character vector of canonical UniProt accessions representing
#'   the full set of detectable proteins (background for the test). Providing
#'   the complete reviewed human proteome gives better-calibrated p-values than
#'   using only the proteins detected in the experiment.
#' @param min_found Integer. Minimum number of annotated complex members that
#'   must be detected to keep a complex (default: \code{2}).
#' @param min_total Integer. Minimum number of annotated members a complex must
#'   have in the lookup table to be tested (default: \code{3}).
#' @param alpha Numeric. Adjusted p-value threshold for retention
#'   (default: \code{0.05}).
#' @param adjust_method Character. Multiple-testing correction method passed to
#'   \code{\link[stats]{p.adjust}} (default: \code{"BH"}).
#'
#' @return List containing:
#'   \itemize{
#'     \item keep: Significantly enriched complexes
#'     \item pruned: Non-significant complexes
#'     \item summary: data.frame with enrichment statistics
#'   }
#'
#' @export
prune_complexes_stat <- function(
  all_complexes,
  lookup_table,
  id_col_data,
  id_col_lookup,
  universe = NULL,
  min_found = 2,
  min_total = 3,
  alpha = 0.05,
  adjust_method = "BH"
) {
  # Deduplicate lookup table
  lookup_table <- lookup_table[
    !duplicated(paste(lookup_table$complex_id, lookup_table[[id_col_lookup]])),
  ]

  # Define universe: all proteins that could be detected
  if (is.null(universe)) {
    # Extract all unique protein IDs from all complexes (handling multiple IDs per row)
    universe <- unique(unlist(lapply(all_complexes, function(df) {
      unlist(lapply(df[[id_col_data]], function(x) {
        if (is.na(x) || x == "") {
          return(character(0))
        }
        trimws(strsplit(x, ";")[[1]])
      }))
    })))
    warning(
      "Universe not provided. Using union of found proteins. For better calibration, provide the full set of measurable proteins."
    )
  }
  universe <- unique(universe[!is.na(universe) & universe != ""])
  N <- length(universe)

  pruned_complexes <- list()
  keep_complexes <- list()
  summary_list <- list()

  for (i in seq_along(all_complexes)) {
    data <- all_complexes[[i]]
    complex_name <- names(all_complexes)[i]

    # Get complex IDs for this complex
    complex_ids <- unique(data$complex_id)

    # Get all annotated proteins for this complex from lookup table (within universe)
    # Normalize: X-1 (isoform-1) is treated as canonical X
    annotated <- unique(sub("-1$", "", lookup_table[[id_col_lookup]][
      lookup_table$complex_name == complex_name &
        !is.na(lookup_table[[id_col_lookup]])
    ]))
    annotated <- intersect(annotated, universe)
    K <- length(annotated)

    # Get proteins found in data (within universe)
    # Handle multiple IDs per row separated by semicolons; normalize -1 isoforms
    found <- unique(sub("-1$", "", unlist(lapply(data[[id_col_data]], function(x) {
      if (is.na(x) || x == "") {
        return(character(0))
      }
      trimws(strsplit(x, ";")[[1]])
    }))))
    found <- intersect(found, universe)
    n <- length(found)

    # Calculate overlap
    k <- length(intersect(found, annotated))

    # Hypergeometric test (right-tail: enrichment)
    p_value <- if (N > 0 && K > 0 && n > 0) {
      phyper(k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
    } else {
      1
    }

    # Count total peptides
    num_peptides <- nrow(data)

    # Store summary info
    summary_list[[i]] <- data.frame(
      complex_name = complex_name,
      complex_ids = paste(complex_ids, collapse = ";"),
      proteins_found = k,
      proteins_total = K,
      query_size = n,
      universe_size = N,
      peptides_total = num_peptides,
      p_value = p_value,
      kept = FALSE,
      stringsAsFactors = FALSE
    )

    # Store in pruned by default
    pruned_complexes[[complex_name]] <- data
  }

  # Combine summary into data frame
  summary_df <- do.call(rbind, summary_list)

  # Adjust p-values for multiple testing
  summary_df$p_adj <- p.adjust(summary_df$p_value, method = adjust_method)

  # Determine which to keep based on criteria
  keep_idx <- summary_df$proteins_found >= min_found &
    summary_df$proteins_total >= min_total &
    summary_df$p_adj <= alpha

  # Move kept complexes from pruned to keep
  for (i in which(keep_idx)) {
    complex_name <- summary_df$complex_name[i]
    keep_complexes[[complex_name]] <- pruned_complexes[[complex_name]]
    pruned_complexes[[complex_name]] <- NULL
    summary_df$kept[i] <- TRUE
  }

  message(paste("Kept", length(keep_complexes), "complexes"))
  message(paste("Pruned", length(pruned_complexes), "complexes"))

  list(
    keep = keep_complexes,
    pruned = pruned_complexes,
    summary = summary_df
  )
}

#' Resolve Ambiguous Peptides Using Complex Membership
#'
#' For peptide rows that map to multiple proteins (semicolon-separated
#' accessions), checks how many of those proteins are annotated members of the
#' complex the row belongs to (via the \code{complex_id} column added by
#' \code{\link{find_complexes}}):
#'
#' \itemize{
#'   \item Exactly one accession in the complex: the row is kept. The accession
#'     column is updated to that single relevant accession, the position column
#'     (if present) is trimmed to the matching \code{"ACC [start-end]"} entry,
#'     and the gene name column (if present) is updated to the corresponding
#'     gene.
#'   \item Multiple accessions in the complex: the peptide is truly ambiguous
#'     within the complex context and the row is removed.
#'   \item No accession in the complex: the row is removed.
#' }
#'
#' Isoform normalization follows the project-wide rule: isoform-1 accessions
#' (e.g. \code{Q13936-1}) are treated as equivalent to the canonical base
#' accession (\code{Q13936}).
#'
#' @param data data.table of complex peptide rows containing a \code{complex_id}
#'   column (added by \code{\link{find_complexes}}).
#' @param accession_col Character. Column name for protein accessions
#'   (default: \code{"Accession"}).
#' @param position_col Character or \code{NULL}. Column containing
#'   semicolon-separated position strings of the form
#'   \code{"ACC [start-end]; ACC2 [start-end]"} (default:
#'   \code{"Position in master protein"}). Set to \code{NULL} to skip.
#' @param gene_col Character or \code{NULL}. Column containing
#'   semicolon-separated gene names mirroring the accession order (default:
#'   \code{"Gene name"}). Set to \code{NULL} to skip.
#' @param lookup_table Complex lookup table (same object passed to
#'   \code{\link{find_complexes}}).
#' @param id_col_lookup Character. Column name for accessions in
#'   \code{lookup_table} (default: \code{"uniprot_id"}).
#'
#' @return A data.table with ambiguous rows either resolved or removed.
#'
#' @export
resolve_ambiguous_peptides <- function(
  data,
  accession_col = "Accession",
  position_col = "Position in master protein",
  gene_col = "Gene name",
  lookup_table,
  id_col_lookup = "uniprot_id"
) {
  if (!accession_col %in% colnames(data)) {
    stop(paste("Column", accession_col, "not found in data"))
  }

  norm <- function(a) sub("-1$", "", trimws(a))

  data <- data.table::copy(data)

  has_position <- !is.null(position_col) && position_col %in% colnames(data)
  has_gene <- !is.null(gene_col) && gene_col %in% colnames(data)

  # Get normalized complex member accessions from lookup
  complex_ids <- unique(data$complex_id)
  complex_members <- norm(unique(lookup_table[[id_col_lookup]][
    lookup_table$complex_id %in% complex_ids &
      !is.na(lookup_table[[id_col_lookup]])
  ]))

  is_ambiguous <- grepl(";", data[[accession_col]], fixed = TRUE)

  if (!any(is_ambiguous)) {
    return(data)
  }

  keep <- rep(TRUE, nrow(data))

  for (i in which(is_ambiguous)) {
    accs <- trimws(strsplit(data[[accession_col]][i], ";")[[1]])
    in_complex <- norm(accs) %in% complex_members

    if (sum(in_complex) == 1L) {
      idx <- which(in_complex)
      resolved_acc <- accs[idx]

      data.table::set(data, i, accession_col, resolved_acc)

      # Keep only the position entry for the resolved accession.
      # Each entry has the form "ACC [start-end]", so match on "ACC [".
      if (has_position) {
        pos_str <- data[[position_col]][i]
        if (!is.na(pos_str) && nzchar(pos_str)) {
          pos_parts <- trimws(strsplit(pos_str, ";")[[1]])
          match_prefix <- paste0(resolved_acc, " [")
          matched_pos <- pos_parts[startsWith(pos_parts, match_prefix)]
          if (length(matched_pos) >= 1L) {
            data.table::set(data, i, position_col, matched_pos[1L])
          }
        }
      }

      # Keep only the gene name at the same index as the resolved accession.
      if (has_gene) {
        gene_str <- data[[gene_col]][i]
        if (!is.na(gene_str) && nzchar(gene_str)) {
          gene_parts <- trimws(strsplit(gene_str, ";")[[1]])
          if (idx <= length(gene_parts)) {
            data.table::set(data, i, gene_col, gene_parts[idx])
          }
        }
      }
    } else {
      keep[i] <- FALSE
    }
  }

  data[keep]
}
