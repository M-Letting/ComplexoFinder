# Initialize time
start_time <- Sys.time()
# Set to not show warnings
options(warn = -1)

# Load libraries
library(data.table, warn.conflicts = FALSE)

# Set working directory to Q-Benchmark
setwd(here::here("Q-Benchmark"))

################################################################################
### Create ProteoMaker datasets ################################################
################################################################################

# Specify datasets
datasets <- list(
  "LowNA" = "./02-DCF/00-Data/01-Unprepared/ProteoMaker_LowNA.csv",
  "LowNoNA" = "./02-DCF/00-Data/01-Unprepared/ProteoMaker_LowNoNA.csv",
  "MedNA" = "./02-DCF/00-Data/01-Unprepared/ProteoMaker_MedNA.csv",
  "MedNoNA" = "./02-DCF/00-Data/01-Unprepared/ProteoMaker_MedNoNA.csv",
  "HighNA" = "./02-DCF/00-Data/01-Unprepared/ProteoMaker_HighNA.csv",
  "HighNoNA" = "./02-DCF/00-Data/01-Unprepared/ProteoMaker_HighNoNA.csv"
)

# Check that files exist and read them
PM_data <- list()
for (name in names(datasets)) {
  file_path <- datasets[[name]]
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  PM_data[[name]] <- fread(file_path)
  cat(
    "Data read from",
    file_path,
    "with",
    nrow(PM_data[[name]]),
    "rows and",
    ncol(PM_data[[name]]),
    "columns.\n"
  )
}

# Get intensity columns (C_1_R_1 etc.)
intensity_cols <- grep("^C_\\d+_R_\\d+$", names(PM_data[[1]]), value = TRUE)

# Specify columns to keep
columns_to_keep <- c(
  "Peptide",
  "Accession",
  "Start",
  "Stop",
  "Proteoform_ID",
  "Peptidoform",
  intensity_cols
)

# Prepare, clean, and process datasets
for (dataset_name in names(PM_data)) {
  data <- PM_data[[dataset_name]]

  # Intensity columns for this dataset (do not reuse mutated names across loops)
  dataset_intensity_cols <- grep("^C_\\d+_R_\\d+$", names(data), value = TRUE)
  dataset_columns_to_keep <- c(
    "Peptide",
    "Accession",
    "Start",
    "Stop",
    "Proteoform_ID",
    "Peptidoform",
    dataset_intensity_cols
  )

  # Subset columns
  data <- data[, ..dataset_columns_to_keep]

  check_cols <- c("Peptide", "Accession", "Start", "Stop")

  # Remove duplicated values inside "|" separated strings
  clean_pipe_values <- function(x) {
    split_vals <- strsplit(x, "\\|", fixed = FALSE)

    sapply(
      seq_along(split_vals),
      function(i) {
        vals <- split_vals[[i]]

        if (length(vals) == 0 || is.na(x[i]) || x[i] == "") {
          return(x[i])
        }

        vals <- unique(trimws(vals))

        paste(vals, collapse = "|")
      },
      USE.NAMES = FALSE
    )
  }

  data[, (check_cols) := lapply(.SD, clean_pipe_values), .SDcols = check_cols]

  # ----- Create Position column -----

  start_split <- strsplit(data$Start, "\\|")
  stop_split <- strsplit(data$Stop, "\\|")

  data[,
    Position := mapply(
      function(s, e) {
        if (length(s) != length(e)) {
          stop("Mismatched number of Start and Stop values")
        }

        paste(paste(s, e, sep = "-"), collapse = "|")
      },
      start_split,
      stop_split,
      USE.NAMES = FALSE
    )
  ]

  # Remove Start/Stop
  data[, c("Start", "Stop") := NULL]

  setcolorder(
    data,
    c(
      "Peptide",
      "Accession",
      "Position",
      "Proteoform_ID",
      "Peptidoform",
      dataset_intensity_cols
    )
  )

  # Rename Intensity columns fromn C_1_R_1 to C1_R1 etc.
  new_intensity_cols <- gsub(
    "([A-Za-z])_([0-9])",
    "\\1\\2",
    dataset_intensity_cols
  )
  setnames(data, dataset_intensity_cols, new_intensity_cols)

  # Median normalize intensity columns
  for (col in new_intensity_cols) {
    col_median <- median(data[[col]], na.rm = TRUE)
    data[[col]] <- data[[col]] - col_median
  }

  # Update dataset
  PM_data[[dataset_name]] <- data
}

# Save processed datasets
for (dataset_name in names(PM_data)) {
  output_path <- paste0(
    "./02-DCF/00-Data/02-Prepared/ProteoMaker",
    dataset_name,
    "_wide_processed.csv"
  )
  fwrite(PM_data[[dataset_name]], output_path)
}
