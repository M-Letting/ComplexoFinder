library(data.table)

# canonical_accession() is defined in 00-LookupTables/GetLookupTables.R

#' Split a semicolon-separated accession string into its accessions
#'
#' @param accession_string Character. One accession, or several separated by
#'   semicolons.
#' @return Character vector of accessions; empty for \code{NA} or \code{""}.
split_accessions <- function(accession_string) {
  if (is.na(accession_string) || accession_string == "") {
    return(character(0))
  }
  trimws(strsplit(accession_string, ";", fixed = TRUE)[[1]])
}

#' Member position of each lookup row
#'
#' @param lookup_table Complex lookup table.
#' @param lookup_id_column Character. Column in \code{lookup_table} with
#'   canonical accessions.
#'
#' @return Character vector, one per lookup row: the \code{member_slot} column
#'   of the lookup, or the canonical accession when the lookup has no such
#'   column (every protein is then its own position).
lookup_member_slot <- function(lookup_table, lookup_id_column) {
  if ("member_slot" %in% names(lookup_table)) {
    lookup_table$member_slot
  } else {
    canonical_accession(lookup_table[[lookup_id_column]])
  }
}

#' Find Protein Complexes in Dataset
#'
#' Identifies protein complexes in a list of datasets by matching protein
#' accessions against a lookup table (EBI Complex Portal, see
#' \code{build_ebi_lookup()}).
#'
#' Matching is on canonical accessions: isoform and processed-chain suffixes
#' are removed on both sides (\code{canonical_accession()}), so a peptide
#' reported for \code{Q13936-2} belongs to every complex that lists
#' \code{Q13936} in any form. Rows that name several proteins
#' (semicolon-separated) match when any of them is a member.
#'
#' @param data_list Named list of data.tables containing proteomics data, one
#'   per datatype.
#' @param lookup_table data.table with complex annotations: \code{complex_id},
#'   \code{complex_name} and the accession column \code{lookup_id_column}.
#' @param data_id_column Character. Column in the data with protein accessions.
#' @param lookup_id_column Character. Column in \code{lookup_table} with
#'   canonical accessions.
#'
#' @return data.table with one row per complex and matching data row: the
#'   columns \code{complex_id} and \code{complex_name}, followed by the data
#'   columns. A data row appears once for every complex it belongs to.
#'
#' @export
find_complexes <- function(
  data_list,
  lookup_table,
  data_id_column = "Accession",
  lookup_id_column = "uniprot_id"
) {
  complex_members <- unique(data.table(
    complex_id = lookup_table$complex_id,
    complex_name = lookup_table$complex_name,
    .match_id = canonical_accession(lookup_table[[lookup_id_column]])
  ))

  matches_by_datatype <- lapply(data_list, function(datatype_table) {
    datatype_table <- as.data.table(datatype_table)

    # Split multi-accession rows into one row per accession for matching
    row_accessions <- datatype_table[,
      .(.match_id = canonical_accession(split_accessions(get(data_id_column)))),
      by = .(.row_index = seq_len(nrow(datatype_table)))
    ]

    matched_rows <- unique(complex_members[
      row_accessions,
      on = ".match_id",
      nomatch = NULL,
      allow.cartesian = TRUE
    ][, .(complex_id, complex_name, .row_index)])

    cbind(
      matched_rows[, .(complex_id, complex_name)],
      datatype_table[matched_rows$.row_index]
    )
  })

  complex_data <- rbindlist(matches_by_datatype, use.names = TRUE, fill = TRUE)
  setorder(complex_data, complex_name)
  complex_data[]
}

#' Prune Complexes by Detected Members
#'
#' Keeps the complexes of which enough members were detected to analyse them
#' as an assembly:
#'
#' \itemize{
#'   \item at least \code{min_proteins_found} member proteins detected, and
#'   \item at least \code{min_fraction_found} of the complex's members
#'     detected.
#' }
#'
#' Members are counted as positions in the complex (\code{member_slot} in the
#' lookup): a position that one of several paralogs can fill counts once, and
#' is detected when any of the paralogs is. A member is detected when it has a
#' row in \code{complex_data}.
#'
#' @param complex_data data.table as returned by \code{find_complexes()}.
#' @param lookup_table Complex lookup table (same object passed to
#'   \code{find_complexes()}).
#' @param data_id_column Character. Column in \code{complex_data} with
#'   accessions.
#' @param lookup_id_column Character. Column in \code{lookup_table} with
#'   canonical accessions.
#' @param min_proteins_found Integer. Minimum number of detected member
#'   proteins.
#' @param min_fraction_found Numeric in [0, 1]. Minimum fraction of members
#'   detected.
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{keep}: rows of \code{complex_data} for the complexes kept.
#'     \item \code{summary}: data.table with one row per complex in the lookup
#'       that has any data: \code{members_total}, \code{members_found},
#'       \code{fraction_found}, \code{proteins_found}, \code{peptides_total}
#'       and \code{kept}.
#'   }
#'
#' @export
prune_complexes <- function(
  complex_data,
  lookup_table,
  data_id_column = "Accession",
  lookup_id_column = "uniprot_id",
  min_proteins_found = 2,
  min_fraction_found = 0.75
) {
  complex_members <- unique(data.table(
    complex_id = lookup_table$complex_id,
    member_slot = lookup_member_slot(lookup_table, lookup_id_column),
    .match_id = canonical_accession(lookup_table[[lookup_id_column]])
  ))

  # Detected proteins per complex
  detected_proteins <- complex_data[,
    .(.match_id = unique(canonical_accession(unlist(lapply(
      unique(get(data_id_column)),
      split_accessions
    ))))),
    by = complex_id
  ]
  detected_proteins[, found := TRUE]

  member_detection <- detected_proteins[
    complex_members,
    on = c("complex_id", ".match_id")
  ]
  member_detection[is.na(found), found := FALSE]

  complex_summary <- member_detection[,
    .(
      members_total = uniqueN(member_slot),
      members_found = uniqueN(member_slot[found]),
      proteins_found = uniqueN(.match_id[found])
    ),
    by = complex_id
  ]
  complex_summary[, fraction_found := members_found / members_total]

  peptide_counts <- complex_data[,
    .(complex_name = complex_name[1], peptides_total = .N),
    by = complex_id
  ]
  complex_summary <- complex_summary[peptide_counts, on = "complex_id"]
  complex_summary[,
    kept := proteins_found >= min_proteins_found &
      fraction_found >= min_fraction_found
  ]
  setcolorder(complex_summary, c("complex_id", "complex_name"))
  setorder(complex_summary, complex_name)

  message(paste("Kept", sum(complex_summary$kept), "complexes"))
  message(paste("Pruned", sum(!complex_summary$kept), "complexes"))

  list(
    keep = complex_data[
      complex_id %in% complex_summary[kept == TRUE, complex_id]
    ],
    summary = complex_summary[]
  )
}

#' Resolve Ambiguous Peptides Using Complex Membership
#'
#' For peptide rows that map to multiple proteins (semicolon-separated
#' accessions), checks how many members of the complex the row belongs to
#' those proteins are:
#'
#' \itemize{
#'   \item One member: the row is kept. The accession column is set to that
#'     accession, the position column (if present) is trimmed to the matching
#'     \code{"ACC [start-end]"} entry, and the gene name column (if present) to
#'     the corresponding gene.
#'   \item Several members: the peptide is ambiguous within the complex and
#'     the row is removed.
#' }
#'
#' A member is a position in the complex (\code{member_slot} in the lookup).
#' Listed accessions that are forms of one protein (isoforms), or paralogs
#' that are alternatives for the same position (e.g. the three calmodulin
#' genes), are therefore one member, and the first listed is used.
#'
#' @param complex_data data.table as returned by \code{find_complexes()} or
#'   the \code{keep} element of \code{prune_complexes()}.
#' @param lookup_table Complex lookup table.
#' @param accession_column Character. Column name for protein accessions.
#' @param position_column Character or \code{NULL}. Column containing
#'   semicolon-separated position strings of the form
#'   \code{"ACC [start-end]; ACC2 [start-end]"}. Set to \code{NULL} to skip.
#' @param gene_column Character or \code{NULL}. Column containing
#'   semicolon-separated gene names mirroring the accession order. Set to
#'   \code{NULL} to skip.
#' @param lookup_id_column Character. Column with canonical accessions in
#'   \code{lookup_table}.
#'
#' @return data.table with ambiguous rows either resolved or removed.
#'
#' @export
resolve_ambiguous_peptides <- function(
  complex_data,
  lookup_table,
  accession_column = "Accession",
  position_column = "Position in master protein",
  gene_column = "Gene name",
  lookup_id_column = "uniprot_id"
) {
  if (!accession_column %in% colnames(complex_data)) {
    stop(paste("Column", accession_column, "not found in data"))
  }

  resolved_data <- data.table::copy(complex_data)
  has_position <- !is.null(position_column) &&
    position_column %in% colnames(resolved_data)
  has_gene <- !is.null(gene_column) && gene_column %in% colnames(resolved_data)

  is_ambiguous <- grepl(";", resolved_data[[accession_column]], fixed = TRUE)
  if (!any(is_ambiguous)) {
    return(resolved_data)
  }

  # Per complex: member position of each member protein
  member_slots_by_complex <- split(
    setNames(
      lookup_member_slot(lookup_table, lookup_id_column),
      canonical_accession(lookup_table[[lookup_id_column]])
    ),
    lookup_table$complex_id
  )

  keep_row <- rep(TRUE, nrow(resolved_data))

  for (row in which(is_ambiguous)) {
    listed_accessions <- split_accessions(resolved_data[[accession_column]][row])
    member_slots <- member_slots_by_complex[[resolved_data$complex_id[row]]][
      canonical_accession(listed_accessions)
    ]
    is_member <- !is.na(member_slots)

    if (uniqueN(member_slots[is_member]) == 1L) {
      first_member <- which(is_member)[1L]
      resolved_accession <- listed_accessions[first_member]

      data.table::set(resolved_data, row, accession_column, resolved_accession)

      # Keep only the position entry for the resolved accession.
      # Each entry has the form "ACC [start-end]", so match on "ACC [".
      if (has_position) {
        position_string <- resolved_data[[position_column]][row]
        if (!is.na(position_string) && nzchar(position_string)) {
          position_entries <- trimws(
            strsplit(position_string, ";", fixed = TRUE)[[1]]
          )
          matching_entries <- position_entries[startsWith(
            position_entries,
            paste0(resolved_accession, " [")
          )]
          if (length(matching_entries) >= 1L) {
            data.table::set(
              resolved_data, row, position_column, matching_entries[1L]
            )
          }
        }
      }

      # Keep only the gene name at the same index as the resolved accession.
      if (has_gene) {
        gene_string <- resolved_data[[gene_column]][row]
        if (!is.na(gene_string) && nzchar(gene_string)) {
          listed_genes <- trimws(strsplit(gene_string, ";", fixed = TRUE)[[1]])
          if (first_member <= length(listed_genes)) {
            data.table::set(
              resolved_data, row, gene_column, listed_genes[first_member]
            )
          }
        }
      }
    } else {
      keep_row[row] <- FALSE
    }
  }

  resolved_data[keep_row]
}

#' Add Peptide Identifiers
#'
#' Adds the column that identifies a feature throughout the pipeline:
#' accession, position, modifications and datatype, joined by \code{"_"}. Two
#' features can share all four (the modification field does not always tell
#' them apart), so repeats within a complex get a numeric suffix (\code{.1},
#' \code{.2}, ...) and every identifier is unique within its complex.
#'
#' @param complex_data data.table of complex peptide rows.
#' @param id_part_columns Character vector. Columns joined into the
#'   identifier.
#' @param peptide_column Character. Name of the identifier column to add.
#'
#' @return \code{complex_data} with the identifier column added.
#'
#' @export
add_peptide_ids <- function(
  complex_data,
  id_part_columns = c(
    "Accession",
    "Position in master protein",
    "Modifications in master protein",
    "datatype"
  ),
  peptide_column = "Peptide"
) {
  missing_columns <- setdiff(id_part_columns, names(complex_data))
  if (length(missing_columns) > 0) {
    stop(
      "Cannot create peptide identifier column. Missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  identified_data <- data.table::copy(complex_data)
  identified_data[,
    (peptide_column) := do.call(paste, c(.SD, sep = "_")),
    .SDcols = id_part_columns
  ]
  identified_data[,
    (peptide_column) := make.unique(get(peptide_column)),
    by = complex_id
  ]
  identified_data[]
}
