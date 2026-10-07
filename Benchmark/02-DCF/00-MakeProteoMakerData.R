#' #############################################################################
#' Module for generating simulated datasets with ProteoMaker
#' #############################################################################
#'
#' This script generates six simulated datasets using the ProteoMaker package,
#' which are then used for benchmarking in the PeCorA analysis. Each dataset is
#' generated with specific parameters to create varying levels of missing
#' values and proteoform complexity.
#'
#' The resulting datasets are saved as CSV files in the "01-Unprepared" directory.
#' #############################################################################

# Packages
suppressMessages(library(devtools, warn.conflicts = FALSE))
# BiocManager::install("cleaver")
suppressMessages(library(cleaver, warn.conflicts = FALSE))
# devtools::install_github("rr-2/PeptideRanger")
suppressMessages(library(PeptideRanger, warn.conflicts = FALSE))
# devtools::install_github("computproteomics/ProteoMaker")
suppressMessages(library(ProteoMaker, warn.conflicts = FALSE))

# Set working directory to Q-Benchmark
setwd(here::here("Q-Benchmark/"))

# Set output directory
output_dir <- "./02-DiscordantPeptidesProteoforms/00-Data/01-Unprepared"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set seed for reproducibility
set.seed(42)

#################################
## Function to convert columns ##
#################################
convert_cols_to_pipe_df <- function(df, cols, null_string = "NULL") {
  convert_col <- function(col) {
    vapply(
      col,
      function(x) {
        if (is.null(x) || length(x) == 0) {
          return(null_string)
        }
        if (length(x) == 1 && is.na(x)) {
          return(NA_character_)
        }
        paste(as.character(x), collapse = "|")
      },
      character(1)
    )
  }

  out <- df
  for (col in cols) {
    out[[col]] <- convert_col(out[[col]])
  }
  out
}

##################################################################
## Set configurations and parameters for all simulated datasets ##
##################################################################
# ProteoMaker configurations
proteomaker_config <- set_proteomaker(
  fastaFilePath = system.file("Proteomes", package = "ProteoMaker"),
  resultFilePath = paste0(tempdir(), "/SimulatedDataSets"),
  cores = 4,
  clusterType = "PSOCK",
  runStatTests = TRUE,
  calcAllBenchmarks = FALSE
)

Param <- def_param()

# Global parameters
Param$Cores <- proteomaker_config$cores
Param$ClusterType <- proteomaker_config$clusterType

# Experimental design
Param$paramGroundTruth$NumCond <- 10
Param$paramGroundTruth$NumReps <- c(3:5)

# Ground truth parameters
Param$paramGroundTruth$FracModProt <- 0.5
Param$paramGroundTruth$PTMTypes <- list(mods = "ph")
Param$paramGroundTruth$PTMTypesMass <- list(m1 = c(ph = 79.966331))
Param$paramGroundTruth$PTMTypesDistr <- list(m1 = c(ph = 1))
Param$paramGroundTruth$PTMMultipleLambda <- 0.1
Param$paramGroundTruth$ModifiableResidues <- list(
  m1 = list(ph = c("S", "T", "Y"))
)
Param$paramGroundTruth$ModifiableResiduesDistr <- list(
  m1 = list(ph = c(0.86, 0.13, 0.01))
)

# Proteoform Abundance Generation parameters
Param$paramProteoformAb$DiffRegFrac <- 0.5

# MS search parameters
Param$paramMSRun$WrongIDs <- 0

#########################
## 1st dataset - LowNA ##
#########################

# Dataset specific parameters
Param$paramMSRun$PercDetectedVal <- 0.5
Param$paramGroundTruth$PropModPerProt <- 1

# Run the simulations
allBs <- run_sims(Param, proteomaker_config)

# Retrieve results
res <- get_simulation(allBs[[1]]$Param, proteomaker_config)

result <- res$StatsPep

# Remove all columns before "Sequence"
result <- result[, which(colnames(result) == "Sequence"):ncol(result)]

# Paste with "|" as seperator for following columns if value of cell is not NA:
# "Peptide", "Start", "Stop", "MC", "Accession", "Proteoform_ID", "PTMPos", "PTMType"
cols_to_paste <- c(
  "Peptide",
  "Start",
  "Stop",
  "MC",
  "Accession",
  "Proteoform_ID",
  "PTMPos",
  "PTMType"
)

# Remove last 2 columns "min1Reg", "allReg"
result <- result[, -which(colnames(result) %in% c("min1Reg", "allReg"))]

result <- convert_cols_to_pipe_df(result, cols_to_paste)

# Save as CSV file
write.csv(
  result,
  file = "./02-DiscordantPeptidesProteoforms/00-Data/01-Unprepared/ProteoMaker_LowNA.csv",
  row.names = FALSE
)

###########################
## 2nd dataset - LowNoNA ##
###########################
# Dataset specific parameters
Param$paramMSRun$PercDetectedVal <- 1
Param$paramGroundTruth$PropModPerProt <- 1

# Run the simulations
allBs <- run_sims(Param, proteomaker_config)

# Retrieve results
res <- get_simulation(allBs[[1]]$Param, proteomaker_config)

result <- res$StatsPep

# Remove all columns before "Sequence"
result <- result[, which(colnames(result) == "Sequence"):ncol(result)]

# Paste with "|" as seperator for following columns if value of cell is not NA:
# "Peptide", "Start", "Stop", "MC", "Accession", "Proteoform_ID", "PTMPos", "PTMType"
cols_to_paste <- c(
  "Peptide",
  "Start",
  "Stop",
  "MC",
  "Accession",
  "Proteoform_ID",
  "PTMPos",
  "PTMType"
)

# Remove last 2 columns "min1Reg", "allReg"
result <- result[, -which(colnames(result) %in% c("min1Reg", "allReg"))]

result <- convert_cols_to_pipe_df(result, cols_to_paste)

# Save as CSV file
write.csv(
  result,
  file = "./02-DiscordantPeptidesProteoforms/00-Data/01-Unprepared/ProteoMaker_LowNoNA.csv",
  row.names = FALSE
)

#########################
## 3rd dataset - MedNA ##
#########################
# Dataset specific parameters
Param$paramMSRun$PercDetectedVal <- 0.5
Param$paramGroundTruth$PropModPerProt <- 5

# Run the simulations
allBs <- run_sims(Param, proteomaker_config)

# Retrieve results
res <- get_simulation(allBs[[1]]$Param, proteomaker_config)

result <- res$StatsPep

# Remove all columns before "Sequence"
result <- result[, which(colnames(result) == "Sequence"):ncol(result)]

# Paste with "|" as seperator for following columns if value of cell is not NA:
# "Peptide", "Start", "Stop", "MC", "Accession", "Proteoform_ID", "PTMPos", "PTMType"
cols_to_paste <- c(
  "Peptide",
  "Start",
  "Stop",
  "MC",
  "Accession",
  "Proteoform_ID",
  "PTMPos",
  "PTMType"
)

# Remove last 2 columns "min1Reg", "allReg"
result <- result[, -which(colnames(result) %in% c("min1Reg", "allReg"))]

result <- convert_cols_to_pipe_df(result, cols_to_paste)

# Save as CSV file
write.csv(
  result,
  file = "./02-DiscordantPeptidesProteoforms/00-Data/01-Unprepared/ProteoMaker_MedNA.csv",
  row.names = FALSE
)

###########################
## 4th dataset - MedNoNA ##
###########################
# Dataset specific parameters
Param$paramMSRun$PercDetectedVal <- 1
Param$paramGroundTruth$PropModPerProt <- 5

# Run the simulations
allBs <- run_sims(Param, proteomaker_config)

# Retrieve results
res <- get_simulation(allBs[[1]]$Param, proteomaker_config)

result <- res$StatsPep

# Remove all columns before "Sequence"
result <- result[, which(colnames(result) == "Sequence"):ncol(result)]

# Paste with "|" as seperator for following columns if value of cell is not NA:
# "Peptide", "Start", "Stop", "MC", "Accession", "Proteoform_ID", "PTMPos", "PTMType"
cols_to_paste <- c(
  "Peptide",
  "Start",
  "Stop",
  "MC",
  "Accession",
  "Proteoform_ID",
  "PTMPos",
  "PTMType"
)

# Remove last 2 columns "min1Reg", "allReg"
result <- result[, -which(colnames(result) %in% c("min1Reg", "allReg"))]

result <- convert_cols_to_pipe_df(result, cols_to_paste)

# Save as CSV file
write.csv(
  result,
  file = "./02-DiscordantPeptidesProteoforms/00-Data/01-Unprepared/ProteoMaker_MedNoNA.csv",
  row.names = FALSE
)

##########################
## 5th dataset - HighNA ##
##########################
# Dataset specific parameters
Param$paramMSRun$PercDetectedVal <- 0.5
Param$paramGroundTruth$PropModPerProt <- 20

# Run the simulations
allBs <- run_sims(Param, proteomaker_config)

# Retrieve results
res <- get_simulation(allBs[[1]]$Param, proteomaker_config)

result <- res$StatsPep

# Remove all columns before "Sequence"
result <- result[, which(colnames(result) == "Sequence"):ncol(result)]

# Paste with "|" as seperator for following columns if value of cell is not NA:
# "Peptide", "Start", "Stop", "MC", "Accession", "Proteoform_ID", "PTMPos", "PTMType"
cols_to_paste <- c(
  "Peptide",
  "Start",
  "Stop",
  "MC",
  "Accession",
  "Proteoform_ID",
  "PTMPos",
  "PTMType"
)

# Remove last 2 columns "min1Reg", "allReg"
result <- result[, -which(colnames(result) %in% c("min1Reg", "allReg"))]

result <- convert_cols_to_pipe_df(result, cols_to_paste)

# Save as CSV file
write.csv(
  result,
  file = "./02-DiscordantPeptidesProteoforms/00-Data/01-Unprepared/ProteoMaker_HighNA.csv",
  row.names = FALSE
)

############################
## 6th dataset - HighNoNA ##
############################
# Dataset specific parameters
Param$paramMSRun$PercDetectedVal <- 1
Param$paramGroundTruth$PropModPerProt <- 20

# Run the simulations
allBs <- run_sims(Param, proteomaker_config)

# Retrieve results
res <- get_simulation(allBs[[1]]$Param, proteomaker_config)

result <- res$StatsPep

# Remove all columns before "Sequence"
result <- result[, which(colnames(result) == "Sequence"):ncol(result)]

# Paste with "|" as seperator for following columns if value of cell is not NA:
# "Peptide", "Start", "Stop", "MC", "Accession", "Proteoform_ID", "PTMPos", "PTMType"
cols_to_paste <- c(
  "Peptide",
  "Start",
  "Stop",
  "MC",
  "Accession",
  "Proteoform_ID",
  "PTMPos",
  "PTMType"
)

# Remove last 2 columns "min1Reg", "allReg"
result <- result[, -which(colnames(result) %in% c("min1Reg", "allReg"))]

result <- convert_cols_to_pipe_df(result, cols_to_paste)

# Save as CSV file
write.csv(
  result,
  file = "./02-DiscordantPeptidesProteoforms/00-Data/01-Unprepared/ProteoMaker_HighNoNA.csv",
  row.names = FALSE
)
