################################################################################
## Lookup tables                                                              ##
## Builds the EBI Complex Portal lookup used to assign peptides to complexes  ##
################################################################################

# Source this file, then call get_ebi_lookup(). The lookup is built once from
# the Complex Portal table and cached next to this script; delete the cached
# .rds (or call with rebuild = TRUE) to refresh it.
#
# Paths are relative to the project root.
#
# Sources
#   EBI Complex Portal, human complexes:
#     https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv
#   Gene names (reviewed human proteome, accession -> gene name):
#     00-LookupTables/human_proteome.rds, as built by
#     ComplexoFinder/00-LookupTables/01-GetLookupTables.R

library(data.table)

lookup_dir <- "NewImplementation/00-LookupTables"
ebi_complextab_url <- "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv"

#' Canonical UniProt accession
#'
#' Strips isoform suffixes (\code{P12345-2}) and processed-chain suffixes
#' (\code{P12345-PRO_0000012345}), so every form of a protein gets the
#' accession of the protein.
#'
#' @param accessions Character vector of accessions.
#' @return Character vector of canonical accessions.
#' @export
canonical_accession <- function(accessions) {
  sub("-(PRO_)?[0-9]+$", "", trimws(accessions))
}

#' Build the EBI Complex Portal lookup
#'
#' Reads the Complex Portal table and returns one row per complex and member
#' protein.
#'
#' Members are taken from the "Expanded participant list" column, in which
#' complexes that are members of other complexes are replaced by their
#' proteins, and which holds proteins only (no small molecules or RNA).
#'
#' A participant can be a set of alternatives, written \code{[P1,P2,P3](n)}:
#' one position in the complex that any of the listed paralogs can fill. Every
#' alternative becomes a row, and they share one \code{member_slot}, so a set
#' counts as a single member when the detected fraction of a complex is
#' computed (see \code{prune_complexes()}).
#'
#' @param complex_table_file Path to the Complex Portal table (9606.tsv).
#' @param gene_names data.frame with columns \code{Accession} and
#'   \code{Gene_name}, or \code{NULL} to skip gene names.
#'
#' @return data.table with columns
#'   \itemize{
#'     \item \code{complex_id}, \code{complex_name}: \code{complex_name} is
#'       unique; a name shared by several complexes gets its id appended.
#'     \item \code{uniprot_id}: canonical accession of the member.
#'     \item \code{source_id}: the member as written in the Complex Portal
#'       (isoform or processed chain where given).
#'     \item \code{member_slot}: the position this member fills; alternatives
#'       of one set share a slot.
#'     \item \code{Stoichiometry}, \code{gene_name}.
#'   }
#' @export
build_ebi_lookup <- function(complex_table_file, gene_names = NULL) {
  complex_table <- fread(complex_table_file, sep = "\t", header = TRUE, quote = "")
  complexes <- complex_table[, .(
    complex_id = `#Complex ac`,
    complex_name = `Recommended name`,
    participants = `Expanded participant list`
  )]

  # A few recommended names are shared by more than one complex
  complexes[, complexes_with_name := .N, by = complex_name]
  complexes[
    complexes_with_name > 1,
    complex_name := paste0(complex_name, " [", complex_id, "]")
  ]
  complexes[, complexes_with_name := NULL]

  # One row per participant
  members <- complexes[,
    .(participant = unlist(strsplit(participants, "|", fixed = TRUE))),
    by = .(complex_id, complex_name)
  ]
  members[, Stoichiometry := as.integer(sub(".*\\((\\d+)\\)$", "\\1", participant))]
  members[, participant := sub("\\(\\d+\\)$", "", participant)]

  # Sets of alternatives: [P1,P2] -> one row per alternative, one shared slot
  members[, participant := gsub("\\[|\\]", "", participant)]
  members <- members[,
    .(source_id = trimws(unlist(strsplit(participant, ",", fixed = TRUE)))),
    by = .(complex_id, complex_name, participant, Stoichiometry)
  ]
  members[, uniprot_id := canonical_accession(source_id)]

  # The slot is named by the canonical accessions that can fill it, so isoforms
  # of one protein listed separately collapse into one slot
  members[,
    member_slot := paste(sort(unique(uniprot_id)), collapse = ","),
    by = .(complex_id, participant)
  ]

  members <- members[!duplicated(paste(complex_id, uniprot_id))]

  if (!is.null(gene_names)) {
    members[,
      gene_name := gene_names$Gene_name[match(uniprot_id, gene_names$Accession)]
    ]
  } else {
    members[, gene_name := NA_character_]
  }

  members[, .(
    complex_id,
    complex_name,
    uniprot_id,
    source_id,
    member_slot,
    Stoichiometry,
    gene_name
  )]
}

#' Load the EBI Complex Portal lookup, building it if needed
#'
#' Returns the cached lookup if there is one. Otherwise downloads the Complex
#' Portal table (unless it is already on disk), builds the lookup and caches
#' it.
#'
#' @param rebuild Logical. Download the Complex Portal table again and rebuild
#'   the lookup even if a cached copy exists.
#'
#' @return data.table, see \code{build_ebi_lookup()}.
#' @export
get_ebi_lookup <- function(rebuild = FALSE) {
  lookup_file <- file.path(lookup_dir, "ebi_cp_lookup.rds")
  complex_table_file <- file.path(lookup_dir, "9606_complextab.tsv")

  if (file.exists(lookup_file) && !rebuild) {
    return(readRDS(lookup_file))
  }

  if (!file.exists(complex_table_file) || rebuild) {
    cat("Downloading EBI Complex Portal table...\n")
    download.file(ebi_complextab_url, complex_table_file, quiet = TRUE)
  }

  gene_names_file <- file.path(lookup_dir, "human_proteome.rds")
  gene_names <- if (file.exists(gene_names_file)) readRDS(gene_names_file) else NULL

  lookup <- build_ebi_lookup(complex_table_file, gene_names)
  saveRDS(lookup, lookup_file)
  cat(
    "Built EBI Complex Portal lookup:",
    uniqueN(lookup$complex_id),
    "complexes,",
    nrow(lookup),
    "member rows\n"
  )
  lookup
}

#' Protein sequences for the sequence plots
#'
#' Reads the UniProt sequence library.
#'
#' @param sequence_file Path to the sequence library (.rds with columns
#'   \code{accession} and \code{sequence}).
#'
#' @return Named character vector (accession -> sequence), or \code{NULL} if
#'   the file does not exist.
#' @export
get_sequence_lookup <- function(
  sequence_file = "ComplexoFinder/00-LookupTables/all_human_sequences.rds"
) {
  if (!file.exists(sequence_file)) {
    return(NULL)
  }
  sequence_library <- readRDS(sequence_file)
  setNames(sequence_library$sequence, sequence_library$accession)
}
