################################################################################
## Data Loader and Preprocessor                                               ##
## Loads raw Excel files, preprocesses data, and saves as csv files           ##
################################################################################

# Description:
# 1. Loads lookup tables and raw data from Excel files in data_excel/ directory
# 2. Processes each dataset
#  - Translates accessions to gene names
#  - Finds peptide positions in protein sequences
#  - Standardizes column names across datasets
#  - Standardizes column order across datasets
# 3. Saves processed datasets as csv files in ComplexoFinder/00-Data/

# Libraries
library(data.table)
library(readxl)

# Set working directory to project root
setwd(here::here())

cat("Loading and processing Excel data files...\n\n")

# Create output directories if they don't exist
if (!dir.exists("ComplexoFinder/00-Data")) {
  dir.create("ComplexoFinder/00-Data", recursive = TRUE)
}

################################################################################
## 1. Source 01-GetLookupTables and define functions ###########################
################################################################################

try(
  source("ComplexoFinder/00-LookupTables/01-GetLookupTables.R"),
  silent = TRUE
)

# Function to translate accessions to gene names
translate_accessions <- function(accession_vector) {
  sapply(
    accession_vector,
    function(accession_string) {
      if (is.na(accession_string) || accession_string == "") {
        return(NA_character_)
      }

      # Split by semicolon and space
      accessions <- strsplit(accession_string, "; ")[[1]]
      accessions <- trimws(accessions)

      # Map each accession to gene name
      gene_names <- sapply(accessions, function(acc) {
        # Remove isoform suffix if present (e.g., P12345-2 -> P12345)
        base_acc <- sub("-\\d+$", "", acc)
        gene <- gene_lookup[base_acc]
        if (is.null(gene) || is.na(gene) || length(gene) == 0) {
          return(acc) # Keep accession if no mapping found
        }
        return(gene[1])
      })

      # Return concatenated gene names (keep duplicates to maintain mapping)
      paste(gene_names, collapse = "; ")
    },
    USE.NAMES = FALSE
  )
}

# Function to find peptide positions in protein sequences
find_peptide_position <- function(peptide_seq, accession_string) {
  if (
    is.na(peptide_seq) ||
      is.na(accession_string) ||
      peptide_seq == "" ||
      accession_string == ""
  ) {
    return(NA_character_)
  }

  # Clean peptide sequence: remove modifications and non-amino acid characters
  clean_peptide <- gsub("\\[.*?\\]", "", peptide_seq)
  clean_peptide <- gsub("[^A-Z]", "", clean_peptide)

  # Split multiple accessions
  accessions <- strsplit(accession_string, "; ")[[1]]
  accessions <- trimws(accessions)

  # Find position for each accession
  positions <- sapply(accessions, function(acc) {
    # Remove isoform suffix for lookup
    base_acc <- sub("-\\d+$", "", acc)
    protein_seq <- sequence_lookup[base_acc]

    if (is.null(protein_seq) || is.na(protein_seq) || protein_seq == "") {
      return(paste0(acc, " [NOT FOUND]"))
    }

    # Find peptide in protein sequence
    pos <- regexpr(clean_peptide, protein_seq, fixed = TRUE)

    if (pos[1] == -1) {
      return(paste0(acc, " [NOT FOUND]"))
    } else {
      start_pos <- pos[1]
      end_pos <- start_pos + attr(pos, "match.length") - 1
      return(paste0(acc, " [", start_pos, "-", end_pos, "]"))
    }
  })

  return(paste(positions, collapse = "; "))
}

################################################################################
## 2. Define Files to Process and Load #########################################
################################################################################

excel_files <- list(
  NM = "data_excel/CBO_NM_data_no_imputation.xlsx",
  NMpeptide = "data_excel/MRL_SBE_CBO_TMT_NM_peptide_groups.xlsx",
  Phospho = "data_excel/CBO_TMTset1-12_Phospho_peptide_groups.xlsx",
  LysAc = "data_excel/CBO_TMTset1-12_LysAcetyl_peptide_groups.xlsx",
  Deglyco = "data_excel/CBO_TMTset1-12_Deglyco_peptide_groups.xlsx",
  RmCys = "data_excel/CBO_TMTset1-12_RmCys_peptide_groups.xlsx",
  FreeCys = "data_excel/CBO_TMTset1-12_FreeCys_peptide_groups.xlsx"
)

# Check all files exist
missing_files <- names(excel_files)[!file.exists(unlist(excel_files))]
if (length(missing_files) > 0) {
  cat("Missing Excel files for: ", paste(missing_files, collapse = ", "), "\n")
  stop(
    "Please ensure all required Excel files are in the data_excel/ directory."
  )
}

cat("\nLoading Excel files...\n")

all_data <- list()

for (name in names(excel_files)) {
  cat("  Loading", name, "from", basename(excel_files[[name]]), "...\n")

  # Read Excel file
  data <- as.data.table(read_excel(excel_files[[name]]))

  all_data[[name]] <- data
  cat("    Loaded", nrow(data), "rows,", ncol(data), "columns\n")
}

################################################################################
## 3. Process Each Dataset #####################################################
################################################################################

cat("\nProcessing datasets...\n")

processed_data <- list()

################################
## Process NM (Protein-level) ##
################################
cat("\n  Processing NM...\n")

NM <- copy(all_data$NM)

# Identify measurement columns
measurement_cols <- names(NM)[grepl("H9|IMR90", names(NM))]

# Rename columns to standard format: CellType_Batch_Day
# Extract pattern from original column names
hmm <- sub("^[^,]*,\\s*", "", measurement_cols)
hmm_split <- do.call(rbind, strsplit(hmm, ",\\s*"))
hmm_split <- hmm_split[, 2:4]
hmm_split[, 2] <- sub("_.*", "", hmm_split[, 2])
hmm_split[, 3] <- paste0("D", hmm_split[, 3])
hmm_split[, 1] <- paste0("B", hmm_split[, 1])
new_colnames <- paste(hmm_split[, 2], hmm_split[, 1], hmm_split[, 3], sep = "_")
setnames(NM, old = measurement_cols, new = new_colnames)
measurement_cols <- new_colnames

# Filter genes with NA gene names
NM <- NM[!is.na(`Gene name`)]

# Add metadata columns for consistency
NM[, `Position in master protein` := NA_character_]
NM[, `Modifications in master protein` := NA_character_]

# Keep relevant columns
keep_cols <- c(
  "Accession",
  "Gene name",
  "Position in master protein",
  "Modifications in master protein",
  measurement_cols
)
NM <- NM[, ..keep_cols]

# Add datatype column
NM[, datatype := "NM"]

processed_data$NM <- NM
cat("    Final dimensions:", nrow(NM), "rows x", ncol(NM), "columns\n")

#####################################################
## Process NMpeptide (Peptide-level, non-modified) ##
#####################################################
cat("\n  Processing NMpeptide...\n")

NMpeptide <- copy(all_data$NMpeptide)

# Select relevant columns
cols_to_keep <- c(
  "Sequence",
  "Modifications",
  "# PSMs",
  "# Proteins",
  "Master Protein Accessions",
  "Master Protein Descriptions",
  grep("\\(Scaled\\)", colnames(NMpeptide), value = TRUE)
)
NMpeptide <- NMpeptide[, ..cols_to_keep]

# Remove Pool and N/A columns
cols_to_remove <- grep(
  "Pool, Pool|N/A, N/A|n/a, n/a",
  names(NMpeptide),
  value = TRUE
)
if (length(cols_to_remove) > 0) {
  NMpeptide[, (cols_to_remove) := NULL]
}

# Identify measurement columns
measurement_cols <- names(NMpeptide)[grepl("H9|IMR90", names(NMpeptide))]

# Rename columns to standard format
new_col_names <- sub("^[^,]*,\\s*", "", measurement_cols)
new_col_names <- do.call(rbind, strsplit(new_col_names, ",\\s*"))
new_col_names <- new_col_names[, 2:4]
new_col_names[, 2] <- sub("_.*", "", new_col_names[, 2])
new_col_names[, 3] <- paste0("D", new_col_names[, 3])
new_col_names[, 1] <- paste0("B", new_col_names[, 1])
new_col_names <- paste(
  new_col_names[, 2],
  new_col_names[, 1],
  new_col_names[, 3],
  sep = "_"
)
setnames(NMpeptide, old = measurement_cols, new = new_col_names)
measurement_cols <- new_col_names

# Rename accession columns
setnames(
  NMpeptide,
  old = "Master Protein Accessions",
  new = "Accession"
)

# Translate Accessions to Gene names
NMpeptide[, `Gene name` := translate_accessions(Accession)]

# Find peptide positions in protein sequences
cat("    Finding peptide positions in protein sequences...\n")
total_rows <- nrow(NMpeptide)
position_vector <- character(total_rows)

for (i in 1:total_rows) {
  if (i %% 10000 == 0) {
    cat("      Processing row", i, "of", total_rows, "\n")
  }
  position_vector[i] <- find_peptide_position(
    NMpeptide$Sequence[i],
    NMpeptide$Accession[i]
  )
}

NMpeptide[, `Position in master protein` := position_vector]
cat("    Peptide positions found!\n")

# Handle Modifications column
if ("Modifications" %in% names(NMpeptide)) {
  setnames(NMpeptide, "Modifications", "Modifications in master protein")
} else {
  NMpeptide[, `Modifications in master protein` := NA_character_]
}

# Set modifications to NA for all NM peptides
NMpeptide[, `Modifications in master protein` := NA_character_]

# Keep relevant columns
keep_cols <- c(
  "Accession",
  "Gene name",
  "Position in master protein",
  "Modifications in master protein",
  measurement_cols
)
NMpeptide <- NMpeptide[, ..keep_cols]

# Add datatype column
NMpeptide[, datatype := "NMpeptide"]

processed_data$NMpeptide <- NMpeptide
cat(
  "    Final dimensions:",
  nrow(NMpeptide),
  "rows x",
  ncol(NMpeptide),
  "columns\n"
)

##########################
## Process PTM datasets ##
##########################
ptm_names <- c("Phospho", "LysAc", "Deglyco", "RmCys", "FreeCys")

for (name in ptm_names) {
  cat("\n  Processing", name, "...\n")

  data <- copy(all_data[[name]])

  # Select relevant columns
  cols_to_keep <- c(
    "Master Protein Accessions",
    "Master Protein Descriptions",
    "Annotated Sequence",
    "Modifications",
    "# Proteins",
    "# PSMs",
    "Positions in Master Proteins",
    "Modifications in Master Proteins",
    grep("\\(Scaled\\)", colnames(data), value = TRUE)
  )
  data <- data[, ..cols_to_keep]

  # Remove Pool and N/A columns
  cols_to_remove <- grep("Pool, Pool|N/A, N/A", colnames(data), value = TRUE)
  if (length(cols_to_remove) > 0) {
    data[, (cols_to_remove) := NULL]
  }

  # Identify measurement columns
  measurement_cols <- grep("\\(Scaled\\)", colnames(data), value = TRUE)

  # Rename columns to standard format
  new_col_names <- sub("^[^,]*,\\s*", "", measurement_cols)
  new_col_names <- do.call(rbind, strsplit(new_col_names, ",\\s*"))
  new_col_names <- new_col_names[, 2:4]
  new_col_names[, 2] <- sub("_.*", "", new_col_names[, 2])
  new_col_names[, 3] <- paste0("D", new_col_names[, 3])
  new_col_names[, 1] <- paste0("B", new_col_names[, 1])
  new_col_names <- paste(
    new_col_names[, 2],
    new_col_names[, 1],
    new_col_names[, 3],
    sep = "_"
  )
  setnames(data, old = measurement_cols, new = new_col_names)
  measurement_cols <- new_col_names

  # Rename accession columns
  setnames(
    data,
    old = c(
      "Master Protein Accessions",
      "Positions in Master Proteins",
      "Modifications in Master Proteins"
    ),
    new = c(
      "Accession",
      "Position in master protein",
      "Modifications in master protein"
    )
  )

  # Translate Accessions to Gene names
  data[, `Gene name` := translate_accessions(Accession)]

  # Keep relevant columns
  keep_cols <- c(
    "Accession",
    "Gene name",
    "Position in master protein",
    "Modifications in master protein",
    measurement_cols
  )
  data <- data[, ..keep_cols]

  # Add datatype column
  data[, datatype := name]

  processed_data[[name]] <- data
  cat("    Final dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
}

################################################################################
## 4. Reorder Columns for Consistency ##########################################
################################################################################

cat("\nReordering columns...\n")

for (name in names(processed_data)) {
  data <- processed_data[[name]]

  # Define measurement columns
  measurement_pattern <- "^(H9|IMR90)_B\\d+_D\\d+$"
  measurement_cols <- grep(measurement_pattern, colnames(data), value = TRUE)

  # Sort measurement columns by: Cell type, Day, Batch
  ordered_measurement_cols <- measurement_cols[order(
    # Cell type (H9 first, then IMR90)
    sub("^(H9|IMR90)_B(\\d+)_D(\\d+)$", "\\1", measurement_cols),
    # Day number
    as.numeric(sub("^(H9|IMR90)_B(\\d+)_D(\\d+)$", "\\3", measurement_cols)),
    # Batch number
    as.numeric(sub("^(H9|IMR90)_B(\\d+)_D(\\d+)$", "\\2", measurement_cols))
  )]

  # Define metadata column order
  metadata_cols <- c(
    "Accession",
    "Gene name",
    "Position in master protein",
    "Modifications in master protein",
    "datatype"
  )

  # Combine metadata and ordered measurements
  new_col_order <- c(metadata_cols, ordered_measurement_cols)

  # Reorder columns
  data <- data[, ..new_col_order]

  # Update the list
  processed_data[[name]] <- data
}

cat("Column reordering complete!\n")

################################################################################
## Save as CSV files in 00-Data ################################################
################################################################################
dataset_names <- c(
  "NM" = "NonModified_Protein",
  "NMpeptide" = "NonModified_Peptide",
  "Phospho" = "Phospho_Peptide",
  "LysAc" = "AcetylatedLysines_Peptide",
  "Deglyco" = "Deglycosylated_Peptide",
  "RmCys" = "ReversiblyModifiedCysteines_Peptide",
  "FreeCys" = "FreeCysteines_Peptide"
)

for (name in names(processed_data)) {
  data <- processed_data[[name]]
  filename <- paste0("ComplexoFinder/00-Data/", dataset_names[name], ".csv")
  fwrite(data, filename)
  cat("Saved", name, "dataset to", filename, "\n")
}

# ################################################################################
# ## 8. Make dataset for comparative testing #####################################
# ################################################################################
# setwd(here::here())
# source("quantification_complex/02-Identification/IdentifyComplexes.R")

# # Choose which dataset to analyze
# dataset <- "H9" # Options: "H9", "IMR90"

# # Verbose output
# verbose <- TRUE # TRUE for detailed output, FALSE for minimal output

# # Include NM (protein level abundance) in analysis
# include_NM <- FALSE # TRUE to include NM column, FALSE to exclude

# if (dataset == "H9") {
#   data_list <- readRDS("quantification_complex/01-DataProcessing/H9_data.rds")
# } else if (dataset == "IMR90") {
#   data_list <- readRDS(
#     "quantification_complex/01-DataProcessing/IMR90_data.rds"
#   )
# } else {
#   stop("Invalid dataset option. Choose 'H9' or 'IMR90'.")
# }

# # Define intensity column patterns for each dataset (after loading data)
# if (dataset == "H9") {
#   intensity_cols <- grep(
#     "^H9_B[0-9]+_D[0-9]+$",
#     names(data_list[[1]]),
#     value = TRUE
#   )
# } else if (dataset == "IMR90") {
#   intensity_cols <- grep("^IMR90", names(data_list[[1]]), value = TRUE)
# } else {
#   stop("Invalid dataset option. Choose 'H9' or 'IMR90'.")
# }

# ################################################################################
# ## Find complexes using EBI annotated list and prune complexes #################
# ################################################################################

# # Lookup table for EBI complexes
# EBI_LT <- readRDS("quantification_complex/00-LookupTables/ebi_cp_lookup.rds")

# # Find complexes using EBI annotated list
# complexes_EBI <- find_complexes(
#   data_list = data_list,
#   lookup_table = EBI_LT,
#   id_col_data = "Gene name",
#   id_col_lookup = "gene_name"
# )

# human_proteome <- readRDS(
#   "quantification_complex/00-LookupTables/human_proteome.rds"
# )
# universe <- unique(human_proteome$Gene_name)

# # Prune complexes using statistical test
# complexes_EBI <- prune_complexes_stat(
#   all_complexes = complexes_EBI$complexes,
#   lookup_table = EBI_LT,
#   id_col_data = "Gene name",
#   id_col_lookup = "gene_name",
#   universe = universe,
#   min_found = 3,
#   min_total = 3,
#   alpha = 0.01,
#   adjust_method = "BH"
# )

# # Remove features with ambigious peptides in all complexes
# for (complex in names(complexes_EBI$keep)) {
#   data <- complexes_EBI$keep[[complex]]

#   data_clean <- remove_ambiguous_peptides(data, accession_col = "Accession")

#   complexes_EBI$keep[[complex]] <- data_clean
# }

# # Cleanup
# rm(data_list, human_proteome, universe)

# ################################################################################
# ## Select complexes to analyze #################################################
# ################################################################################
# complexes <- list(
#   "B-Wich" = "B-WICH chromatin remodelling complex",
#   "CRD-mediated" = "CRD-mediated mRNA stability complex",
#   "Multiaminoacyl-tRNA" = "Multiaminoacyl-tRNA synthetase complex",
#   "Major-Spliceosomal-B" = "Major Spliceosomal B complex",
#   "26S-proteasome" = "26S proteasome complex",
#   "60S-cytosolic-large-ribosomal-subunit" = "60S cytosolic large ribosomal subunit",
#   "Nuclear-pore-complex" = "Nuclear pore complex",
#   "Eukaryotic-translation-initiation-factor-3-complex" = "Eukaryotic translation initiation factor 3 complex",
#   "Dynein-1-complex-variant-1" = "Dynein-1 complex, variant 1",
#   "Intraflagellar-transport-complex-B" = "Intraflagellar transport complex B",
#   "BSSADCRC" = "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
#   "Dynactin-complex" = "Dynactin complex",
#   "Laminin-213" = "Laminin-213 complex"
# )

# # If not including NM, remove features with datatype == NM from complexes_EBI$keep
# if (!include_NM) {
#   for (complex in names(complexes_EBI$keep)) {
#     data <- complexes_EBI$keep[[complex]]
#     if (!is.null(data)) {
#       data <- data[datatype != "NM"]
#       complexes_EBI$keep[[complex]] <- data
#     }
#   }
# }

# for (complex in names(complexes)) {
#   cat("\n=== Complex:", complex, "===\n")
#   complex_name <- complexes[[complex]]
#   data <- complexes_EBI$keep[[complex_name]]
#   if (is.null(data) || nrow(data) == 0) {
#     cat("No data available for this complex after filtering.\n")
#     next
#   }

#   # Combine metadata into a single identifier for each feature (Peptide)
#   data[,
#     Peptide := paste(
#       `Gene name`,
#       Accession,
#       `Position in master protein`,
#       `Modifications in master protein`,
#       datatype,
#       sep = "_"
#     )
#   ]

#   # Get intensity columns
#   intensity_cols <- grep("^(H9|IMR90)_B\\d+_D\\d+$", names(data), value = TRUE)

#   # Reorder columns and keep only complex_id, Peptide and intensity columns
#   data <- data[, c("complex_id", "Peptide", intensity_cols), with = FALSE]

#   # Remove H9_ prefix from intensity columns for better readability
#   new_intensity_cols <- sub("^(H9|IMR90)_", "", intensity_cols)
#   setnames(data, old = intensity_cols, new = new_intensity_cols)

#   # Save the data for this complex
#   fwrite(
#     data,
#     file = paste0(
#       "Q-Benchmark/02-DCF/00-Data/02-Prepared/",
#       complex,
#       "_WideLogNorm.csv"
#     ),
#     row.names = FALSE
#   )
# }
