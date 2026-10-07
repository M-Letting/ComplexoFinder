# ============================================================================
# BENCHMARK DATA LOADERS
# ============================================================================
# Loaders for the two benchmarks with a known ground truth. Each returns a
# list with
#   data, intensity_columns, assembly_column, peptide_column, condition_regex,
#   datatype_column, truth
# with one row per peptide in `data`.
#
#   SWATH-MS     real DIA data (HEK293), 3 days x 7 replicates, in which some
#                peptides of 1,000 proteins were reduced on one day.
#                Truth: the peptide was perturbed.
#                Input: Benchmark/01-FindDCF/data/prepared/
#
#   ProteoMaker  simulated proteoform data, 10 conditions x 3 replicates.
#                Truth: reconstructed from the simulation metadata. A peptide
#                is discordant when its expected profile differs from the most
#                common expected profile of its protein; peptides of a protein
#                with the same expected profile form a truth group.
#                Input: Benchmark/02-DCF/00-Data/01-Unprepared/ and 02-Prepared/
#
# Paths are relative to the project root.
# ============================================================================

suppressMessages({
  library(data.table)
  library(arrow)
})

swath_datasets <- c("1pep", "2pep", "random", "050pep")
proteomaker_datasets <- c("LowNA", "MedNA", "HighNA", "LowNoNA", "MedNoNA", "HighNoNA")

# ---------------------------------------------------------------------------
# SWATH-MS
# ---------------------------------------------------------------------------

#' Load a SWATH-MS dataset with its ground truth
#'
#' Reads the prepared long table, spreads it to one row per peptide and one
#' column per sample, log2-transforms the intensities and centers every sample
#' on its median.
#'
#' @param dataset_name One of `swath_datasets`.
#'
#' @return List with `data` (one row per peptide), `intensity_columns`,
#'   `assembly_column` (the protein), `peptide_column`, `condition_regex`,
#'   `datatype_column` (`NULL`) and `truth` (data.table with `group_id`,
#'   `peptide_id` and `truth_discordant`: the peptide was perturbed).
load_swath_dataset <- function(dataset_name) {
  input_file <- file.path(
    "Benchmark/01-FindDCF/data/prepared",
    paste0("bench_", dataset_name, "_input.feather")
  )
  if (!file.exists(input_file)) stop("Missing SWATH-MS input: ", input_file)

  long_data <- as.data.table(as.data.frame(arrow::read_feather(input_file)))
  long_data[, sample_id := paste(day, filename, sep = "_")]

  wide_data <- dcast(
    long_data,
    protein_id + peptide_id + n_pep + n_perturbed_peptides +
      perturbed_protein + perturbed_peptide + red_fac ~ sample_id,
    value.var = "intensity",
    fun.aggregate = mean
  )

  intensity_columns <- grep("^day[0-9]+_", names(wide_data), value = TRUE)
  for (column in intensity_columns) {
    wide_data[[column]] <- log2(wide_data[[column]])
  }
  for (column in intensity_columns) {
    wide_data[[column]] <- wide_data[[column]] - median(wide_data[[column]], na.rm = TRUE)
  }

  list(
    data = wide_data,
    intensity_columns = intensity_columns,
    assembly_column = "protein_id",
    peptide_column = "peptide_id",
    condition_regex = "^(day[0-9]+)_.*$",
    datatype_column = NULL,
    truth = unique(wide_data[, .(
      group_id = protein_id,
      peptide_id,
      truth_discordant = as.logical(perturbed_peptide)
    )])
  )
}

# ---------------------------------------------------------------------------
# ProteoMaker
# ---------------------------------------------------------------------------

#' Parse a regulation pattern from the ProteoMaker metadata
#'
#' @param pattern_string Character, e.g. \code{"c(-2, 1.5, -2.5)"}.
#' @return Numeric vector, or \code{NULL} for \code{"NULL"}, \code{"NA"} or
#'   \code{""} (an unregulated proteoform).
.parse_pattern <- function(pattern_string) {
  pattern_string <- trimws(pattern_string)
  if (pattern_string %in% c("NULL", "NA", "")) return(NULL)
  suppressWarnings(as.numeric(
    strsplit(gsub("^c\\(|\\)$", "", pattern_string), ",\\s*")[[1]]
  ))
}

#' Expected condition profile of each peptide, from the ProteoMaker metadata
#'
#' A peptide is carried by one or more proteoforms (Proteoform_ID,
#' "|"-separated), each of which has a Regulation_Pattern (a per-condition
#' shape, or NULL = unregulated = flat) scaled by a Regulation_Amplitude
#' (";"-separated, aligned with Proteoform_ID). The peptide's expected profile
#' is the mean of its proteoforms' amplitude x pattern vectors.
#'
#' @param metadata data.table with the columns `Regulation_Pattern` and
#'   `Regulation_Amplitude`, one row per peptide.
#' @param n_conditions Integer. Number of conditions.
#'
#' @return Numeric matrix, peptides x conditions.
build_expected_profiles <- function(metadata, n_conditions) {
  n_peptides <- nrow(metadata)
  expected_profiles <- matrix(0, n_peptides, n_conditions)
  patterns_by_peptide <- strsplit(metadata$Regulation_Pattern, ";", fixed = TRUE)
  amplitudes_by_peptide <- strsplit(metadata$Regulation_Amplitude, ";", fixed = TRUE)

  for (peptide in seq_len(n_peptides)) {
    patterns <- patterns_by_peptide[[peptide]]
    amplitudes <- amplitudes_by_peptide[[peptide]]
    # One row per proteoform that carries the peptide
    proteoform_profiles <- matrix(0, max(length(patterns), 1L), n_conditions)
    for (proteoform in seq_along(patterns)) {
      pattern <- .parse_pattern(patterns[proteoform])
      if (is.null(pattern) || length(pattern) != n_conditions) next
      amplitude <- suppressWarnings(
        as.numeric(amplitudes[min(proteoform, length(amplitudes))])
      )
      if (!is.finite(amplitude)) amplitude <- 1
      proteoform_profiles[proteoform, ] <- amplitude * pattern
    }
    expected_profiles[peptide, ] <- colMeans(
      proteoform_profiles[seq_along(patterns), , drop = FALSE]
    )
  }
  expected_profiles
}

#' Load a ProteoMaker dataset with its ground truth
#'
#' @param dataset_name One of `proteomaker_datasets`.
#' @param profile_tolerance Expected profiles closer than this (log2) count as
#'   the same profile.
#'
#' @return List with `data` (one row per peptide), `intensity_columns`,
#'   `assembly_column` (the protein), `peptide_column`, `condition_regex`,
#'   `datatype_column`, `expected` (numeric matrix, row-aligned with `data`:
#'   each peptide's expected profile) and `truth` (data.table with `group_id`,
#'   `peptide_id`, `truth_discordant` and `truth_group`, an id shared by the
#'   peptides of a protein that have the same expected profile).
load_proteomaker_dataset <- function(dataset_name, profile_tolerance = 1e-6) {
  metadata_file <- file.path(
    "Benchmark/02-DCF/00-Data/01-Unprepared",
    paste0("ProteoMaker_", dataset_name, ".csv")
  )
  data_file <- file.path(
    "Benchmark/02-DCF/00-Data/02-Prepared",
    paste0("ProteoMaker", dataset_name, "_wide_processed.csv")
  )
  if (!file.exists(metadata_file) || !file.exists(data_file)) {
    stop("Missing ProteoMaker input for dataset: ", dataset_name)
  }

  data <- fread(data_file)
  metadata <- fread(metadata_file)
  # The rows of the two files correspond 1:1 in order
  stopifnot(nrow(metadata) == nrow(data))

  intensity_columns <- grep("^C\\d+_R\\d+$", names(data), value = TRUE)
  condition_regex <- "^(C\\d+)_R\\d+$"
  n_conditions <- length(unique(sub(condition_regex, "\\1", intensity_columns)))

  data[, Peptide_ID := paste(Peptide, Position, Proteoform_ID, Peptidoform, sep = "_")]
  data[, datatype := ifelse(grepl("[ph]", Peptidoform, fixed = TRUE), "Phospho", "NM")]

  expected_profiles <- build_expected_profiles(metadata, n_conditions)
  # One string per distinct expected profile
  profile_signature <- apply(
    round(expected_profiles / profile_tolerance) * profile_tolerance,
    1,
    paste,
    collapse = "_"
  )

  truth <- data.table(
    group_id = data$Accession,
    peptide_id = data$Peptide_ID,
    profile_signature = profile_signature
  )
  # The most common expected profile of a protein is its majority behaviour
  truth[,
    truth_discordant := profile_signature !=
      names(sort(table(profile_signature), decreasing = TRUE))[1],
    by = group_id
  ]
  truth[,
    truth_group := paste(group_id, as.integer(factor(profile_signature)), sep = ":"),
    by = group_id
  ]
  truth[, profile_signature := NULL]

  list(
    data = data,
    intensity_columns = intensity_columns,
    assembly_column = "Accession",
    peptide_column = "Peptide_ID",
    condition_regex = condition_regex,
    datatype_column = "datatype",
    truth = truth,
    expected = expected_profiles
  )
}
