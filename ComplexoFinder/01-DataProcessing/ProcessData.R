################################################################################
## Data Processing                                                            ##
## Loads csv files, splits by cell line, filters, transforms, and normalizes  ##
################################################################################

# Description:
# 1. Loads csv files created by DataLoadPreprocess.R
# 2. Splits data into H9 and IMR90 cell line datasets
# 3. Removes features with >60% missing values within each cell line
# 4. Save untransformed data for benchmarking
# 5. Log-transforms intensity data
#   - Note: NMprotein is not log-transformed as it is already log-transformed
# 6. Median-centers columns (normalization)
# 7. Saves processed data lists as RDS files in ComplexoFinder/00-Data/

# Libraries
library(data.table)

# Set working directory
setwd(here::here("ComplexoFinder/"))

# Check if data files exist in csv format
if (
  !file.exists("00-Data/AcetylatedLysines_Peptide.csv") ||
    !file.exists("00-Data/Deglycosylated_Peptide.csv") ||
    !file.exists("00-Data/FreeCysteines_Peptide.csv") ||
    !file.exists("00-Data/NonModified_Peptide.csv") ||
    !file.exists("00-Data/NonModified_Protein.csv") ||
    !file.exists("00-Data/Phospho_Peptide.csv") ||
    !file.exists("00-Data/ReversiblyModifiedCysteines_Peptide.csv")
) {
  # If not, run data loading and prehandle script
  source("01-DataProcessing/DataLoadPreprocess.R")
}

# Load data into list
processed_data <- list(
  LysAC = fread("00-Data/AcetylatedLysines_Peptide.csv"),
  Deglyco = fread("00-Data/Deglycosylated_Peptide.csv"),
  FreeCys = fread("00-Data/FreeCysteines_Peptide.csv"),
  NMpeptide = fread("00-Data/NonModified_Peptide.csv"),
  NMprotein = fread("00-Data/NonModified_Protein.csv"),
  Phospho = fread("00-Data/Phospho_Peptide.csv"),
  RmCys = fread("00-Data/ReversiblyModifiedCysteines_Peptide.csv")
)

################################################################################
## Split Data by Cell Line #####################################################
################################################################################

cat("\nSplitting data by cell line...\n")

# Split into H9 and IMR90 lists
H9_data <- list()
IMR90_data <- list()

for (name in names(processed_data)) {
  dt <- processed_data[[name]]

  # Create H9 dataset (metadata + H9 measurements)
  h9_cols <- grepl("H9", names(dt))
  metadata_cols <- !grepl("H9|IMR90", names(dt))
  H9_data[[name]] <- dt[, .SD, .SDcols = which(h9_cols | metadata_cols)]

  # Create IMR90 dataset (metadata + IMR90 measurements)
  imr90_cols <- grepl("IMR90", names(dt))
  IMR90_data[[name]] <- dt[, .SD, .SDcols = which(imr90_cols | metadata_cols)]

  cat(
    "  ",
    name,
    ": H9 =",
    sum(h9_cols),
    "samples, IMR90 =",
    sum(imr90_cols),
    "samples\n"
  )
}

################################################################################
## Filter Each Dataset by Cell Line (≥40% present within cell line) ############
################################################################################

th <- 0.6

cat(
  "\nFiltering datasets to keep features with ≥",
  th * 100,
  "% non-missing values within each cell line...\n"
)

# Filter H9 data
cat("\nFiltering H9 datasets...\n")
for (name in names(H9_data)) {
  dt <- H9_data[[name]]

  # Identify measurement columns
  measurement_cols <- grep("^H9_B\\d+_D\\d+$", names(dt), value = TRUE)

  if (length(measurement_cols) == 0) {
    next
  }

  n_before <- nrow(dt)

  # Filter: keep features present in ≥40% of samples
  threshold <- th * length(measurement_cols)
  dt[, n_present := rowSums(!is.na(.SD)), .SDcols = measurement_cols]
  dt <- dt[n_present >= threshold]
  dt[, n_present := NULL]

  n_after <- nrow(dt)

  H9_data[[name]] <- dt
  cat("  ", name, ":", n_before, "->", n_after, "rows\n")
}

# Filter IMR90 data
cat("\nFiltering IMR90 datasets...\n")
for (name in names(IMR90_data)) {
  dt <- IMR90_data[[name]]

  # Identify measurement columns
  measurement_cols <- grep("^IMR90_B\\d+_D\\d+$", names(dt), value = TRUE)

  if (length(measurement_cols) == 0) {
    next
  }

  n_before <- nrow(dt)

  # Filter: keep features present in ≥40% of samples
  threshold <- th * length(measurement_cols)
  dt[, n_present := rowSums(!is.na(.SD)), .SDcols = measurement_cols]
  dt <- dt[n_present >= threshold]
  dt[, n_present := NULL]

  n_after <- nrow(dt)

  IMR90_data[[name]] <- dt
  cat("  ", name, ":", n_before, "->", n_after, "rows\n")
}

cat("\nFiltering complete!\n")

################################################################################
## Save untransformed data for benchmarking ####################################
################################################################################

saveRDS(H9_data, "00-Data/H9_data_NoLogNoNorm.rds")
saveRDS(IMR90_data, "00-Data/IMR90_data_NoLogNoNorm.rds")

################################################################################
## Log-transform + Normalize ###################################################
################################################################################

cat("\nLog-transforming and normalizing datasets...\n")

# Log-transform and normalize H9 data
for (name in names(H9_data)) {
  dt <- H9_data[[name]]

  # Identify measurement columns
  measurement_cols <- grep("^H9_B\\d+_D\\d+$", names(dt), value = TRUE)

  if (length(measurement_cols) == 0) {
    next
  }

  # Do not log-transform NMprotein dataset (already log-transformed)
  if (name != "NMprotein") {
    # Set 0 values to NA before log2 transformation
    dt[,
      (measurement_cols) := lapply(.SD, function(x) ifelse(x == 0, NA, x)),
      .SDcols = measurement_cols
    ]

    # Log tranform
    dt[, (measurement_cols) := lapply(.SD, log2), .SDcols = measurement_cols]
  }

  # Normalize by median centering
  dt[,
    (measurement_cols) := lapply(.SD, function(x) x - median(x, na.rm = TRUE)),
    .SDcols = measurement_cols
  ]

  H9_data[[name]] <- dt
  cat("  H9 ", name, ": log2-transformed and median-centered\n")
}

# Log-transform and normalize IMR90 data
for (name in names(IMR90_data)) {
  dt <- IMR90_data[[name]]

  # Identify measurement columns
  measurement_cols <- grep("^IMR90_B\\d+_D\\d+$", names(dt), value = TRUE)

  if (length(measurement_cols) == 0) {
    next
  }

  # Do not log-transform NMprotein dataset (already log-transformed)
  if (name != "NMprotein") {
    # Set 0 values to NA before log transformation
    dt[,
      (measurement_cols) := lapply(.SD, function(x) ifelse(x == 0, NA, x)),
      .SDcols = measurement_cols
    ]

    # Log tranform
    dt[, (measurement_cols) := lapply(.SD, log2), .SDcols = measurement_cols]
  }

  # Normalize by median centering
  dt[,
    (measurement_cols) := lapply(.SD, function(x) x - median(x, na.rm = TRUE)),
    .SDcols = measurement_cols
  ]

  IMR90_data[[name]] <- dt
  cat("  IMR90 ", name, ": log2-transformed and median-centered\n")
}

################################################################################
## Save Data Lists as RDS files ################################################
################################################################################

saveRDS(H9_data, "00-Data/H9_data.rds")
saveRDS(IMR90_data, "00-Data/IMR90_data.rds")

cat("\nSaved data lists:\n")
cat("  - 00-Data/H9_data.rds\n")
cat("  - 00-Data/IMR90_data.rds\n")

cat("\n=== Data loading and processing complete! ===\n")
cat("\nSummary:\n")
cat("  Datasets processed:", length(processed_data), "\n")
cat("  H9 datasets:", length(H9_data), "\n")
cat("  IMR90 datasets:", length(IMR90_data), "\n")
cat("  Output files saved in 00-Data/\n")
