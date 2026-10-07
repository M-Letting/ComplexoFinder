################################################################################
## Data Processing                                                            ##
## Loads csv files, splits by cell line, transforms, and normalizes           ##
################################################################################

# Description:
# 1. Loads the csv files created by ComplexoFinder/01-DataProcessing/
#    DataLoadPreprocess.R (one per datatype, both cell lines)
# 2. Splits data into H9 and IMR90 cell line datasets
# 3. Sets 0 intensities to NA and removes features with no value at all within
#    the cell line
# 4. Saves untransformed data for benchmarking
# 5. Log-transforms intensity data
#   - Note: NMprotein is not log-transformed as it is already log-transformed
# 6. Median-centers columns (normalization)
# 7. Saves processed data lists as RDS files in NewImplementation/00-Data/
#
# No completeness filter is applied here: discover_complexoforms() decides
# which peptides are complete enough to test.
#
# Run from the project root:
#   Rscript NewImplementation/01-DataProcessing/ProcessData.R

# Libraries
library(data.table)

# Set working directory
setwd(here::here())

csv_directory <- "ComplexoFinder/00-Data"
output_directory <- "NewImplementation/00-Data"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

csv_files <- c(
  LysAc = "AcetylatedLysines_Peptide.csv",
  Deglyco = "Deglycosylated_Peptide.csv",
  FreeCys = "FreeCysteines_Peptide.csv",
  NMpeptide = "NonModified_Peptide.csv",
  NMprotein = "NonModified_Protein.csv",
  Phospho = "Phospho_Peptide.csv",
  RmCys = "ReversiblyModifiedCysteines_Peptide.csv"
)

missing_csv_files <- csv_files[!file.exists(file.path(csv_directory, csv_files))]
if (length(missing_csv_files) > 0) {
  stop(
    "Missing csv files in ",
    csv_directory,
    ": ",
    paste(missing_csv_files, collapse = ", "),
    "\nRun ComplexoFinder/01-DataProcessing/DataLoadPreprocess.R first."
  )
}

# Load data into list, one table per datatype
datatype_tables <- lapply(file.path(csv_directory, csv_files), fread)
names(datatype_tables) <- names(csv_files)

################################################################################
## Process one cell line #######################################################
################################################################################

#' Split out one cell line and remove features without any value
#'
#' Keeps the metadata columns and the measurement columns of one cell line,
#' sets intensities of 0 to \code{NA} (a 0 is a missing measurement) and drops
#' the features that have no value in the cell line.
#'
#' @param datatype_table data.table of one datatype, both cell lines.
#' @param cell_line Character. "H9" or "IMR90".
#' @param is_log_transformed Logical. \code{TRUE} if the intensities are
#'   already log-transformed, in which case 0 is a valid value and is kept.
#'
#' @return data.table with metadata and the cell line's measurement columns.
split_cell_line <- function(datatype_table, cell_line, is_log_transformed) {
  measurement_columns <- grep(
    paste0("^", cell_line, "_B\\d+_D\\d+$"),
    names(datatype_table),
    value = TRUE
  )
  metadata_columns <- names(datatype_table)[
    !grepl("^(H9|IMR90)_B\\d+_D\\d+$", names(datatype_table))
  ]
  cell_line_table <- datatype_table[,
    c(metadata_columns, measurement_columns),
    with = FALSE
  ]

  # Set 0 values to NA: a 0 intensity is a missing measurement
  if (!is_log_transformed) {
    cell_line_table[,
      (measurement_columns) := lapply(
        .SD,
        function(intensity) ifelse(intensity == 0, NA, intensity)
      ),
      .SDcols = measurement_columns
    ]
  }

  # Remove features never measured in this cell line
  is_measured <- rowSums(
    !is.na(as.matrix(cell_line_table[, ..measurement_columns]))
  ) > 0
  cell_line_table[is_measured]
}

#' Log2-transform and median-center the measurement columns
#'
#' @param cell_line_table data.table from \code{split_cell_line()}.
#' @param cell_line Character. "H9" or "IMR90".
#' @param is_log_transformed Logical. \code{TRUE} if the intensities are
#'   already log-transformed; they are then only median-centered.
#'
#' @return data.table with every measurement column log2-transformed and
#'   centered on its median.
log_normalize <- function(cell_line_table, cell_line, is_log_transformed) {
  normalized_table <- copy(cell_line_table)
  measurement_columns <- grep(
    paste0("^", cell_line, "_B\\d+_D\\d+$"),
    names(normalized_table),
    value = TRUE
  )

  # Log transform
  if (!is_log_transformed) {
    normalized_table[,
      (measurement_columns) := lapply(.SD, log2),
      .SDcols = measurement_columns
    ]
  }

  # Normalize by median centering
  normalized_table[,
    (measurement_columns) := lapply(
      .SD,
      function(intensity) intensity - median(intensity, na.rm = TRUE)
    ),
    .SDcols = measurement_columns
  ]
  normalized_table
}

for (cell_line in c("H9", "IMR90")) {
  cat("\nProcessing", cell_line, "datasets...\n")

  untransformed_tables <- list()
  normalized_tables <- list()

  for (datatype in names(datatype_tables)) {
    # NMprotein is already log-transformed
    is_log_transformed <- datatype == "NMprotein"

    rows_before <- nrow(datatype_tables[[datatype]])
    untransformed_tables[[datatype]] <- split_cell_line(
      datatype_tables[[datatype]],
      cell_line,
      is_log_transformed
    )
    normalized_tables[[datatype]] <- log_normalize(
      untransformed_tables[[datatype]],
      cell_line,
      is_log_transformed
    )

    cat(
      "  ", datatype, ":", rows_before, "->",
      nrow(untransformed_tables[[datatype]]), "rows\n"
    )
  }

  # Untransformed data for benchmarking
  saveRDS(
    untransformed_tables,
    file.path(output_directory, paste0(cell_line, "_data_NoLogNoNorm.rds"))
  )
  saveRDS(
    normalized_tables,
    file.path(output_directory, paste0(cell_line, "_data.rds"))
  )
}

cat("\nSaved data lists in", output_directory, "\n")
