# Set working directory to folder
setwd(here::here("ComplexoFinder"))

# Check if data files exist, if not run data loading and handling script
if (
  !file.exists("00-Data/H9_data.rds") ||
    !file.exists("00-Data/IMR90_data.rds")
) {
  source("01-DataProcessing/ProcessData.R")
}

# Check if required lookup tables exist, if not run script to create them
if (
  !file.exists("00-LookupTables/ebi_cp_lookup.rds") ||
    !file.exists("00-LookupTables/human_proteome.rds")
) {
  source("00-LookupTables/01-GetLookupTables.R")
}

# Check if output directory exists, if not create it
if (!dir.exists("10-Results")) {
  dir.create("10-Results")
}

# Source functions
source("02-Identification/AssignComplexes.R")
source("02-Identification/IdentifyComplexes.R")
source("03-Constraints/Constraints.R")
source("04-Clustering/VsClust.R")
source("04-Clustering/VsClustWrappers.R")
source("04-Clustering/ClusterSummary.R")
source("05-Quantification/Quantification.R")
source("09-Visualization/ExpressionMembershipPlot.R")
source("09-Visualization/SequencePlot.R")
source("09-Visualization/QuantificationPlot.R")
source("09-Visualization/GroupingProteinPCA.R")

################################################################################
## Options #####################################################################
################################################################################

# Choose which dataset to analyze
dataset <- "H9" # Options: "H9", "IMR90"

# Verbose output
verbose <- TRUE # TRUE for detailed output, FALSE for minimal output

# Include NM (protein level abundance) in analysis
include_NM <- FALSE # TRUE to include protein-level data, FALSE to exclude

# Include constrainsts in clustering
include_constraints <- TRUE # TRUE to include constraints, FALSE to exclude

# Specify complexes to analyze - for all options, see EBI complex lookup table
complexes <- c(
  # "B-WICH chromatin remodelling complex",
  # "CRD-mediated mRNA stability complex",
  # "Multiaminoacyl-tRNA synthetase complex",
  # "Major Spliceosomal B complex",
  # "26S proteasome complex",
  # "60S cytosolic large ribosomal subunit",
  # "Nuclear pore complex",
  "Eukaryotic translation initiation factor 3 complex" #,
  # "Dynein-1 complex, variant 1",
  # "Intraflagellar transport complex B",
  # "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  # "Neuron-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  # "Neural progenitor-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  # "Calcineurin-Calmodulin-AKAP5 complex, gamma-R1 variant",
  # "SNARE complex STX4-SNAP29-SEC22b",
  # "Cortical microtubule stabilization complex, KANK1 variant",
  # "LIFT actin modulation complex",
  # "Dynein-1 complex, variant 4",
  # "WASH complex, variant WASHC1/WASHC2C",
  # "Cav1.2 voltage-gated calcium channel complex, CACNA2D1-CACNB3 variant",
  # "Dynactin complex",
  # "Laminin-213 complex"
)

################################################################################
## Run ComplexoFinder ##########################################################
################################################################################

#############
# Load data #
#############

if (dataset == "H9") {
  data_list <- readRDS("00-Data/H9_data.rds")
} else if (dataset == "IMR90") {
  data_list <- readRDS("00-Data/IMR90_data.rds")
} else {
  stop("Invalid dataset choice. Please choose 'H9' or 'IMR90'.")
}

# Define intensity column patterns for each dataset (after loading data)
if (dataset == "H9") {
  intensity_cols <- grep(
    "^H9_B[0-9]+_D[0-9]+$",
    names(data_list[[1]]),
    value = TRUE
  )
} else if (dataset == "IMR90") {
  intensity_cols <- grep("^IMR90", names(data_list[[1]]), value = TRUE)
} else {
  stop("Invalid dataset option. Choose 'H9' or 'IMR90'.")
}

n_rep <- 3
n_cond <- length(intensity_cols) / n_rep
if (n_cond != as.integer(n_cond)) {
  stop("Cannot infer n_cond: length(intensity_cols) is not divisible by n_rep")
}
n_cond <- as.integer(n_cond)

cond_regex <- if (dataset == "H9") {
  "^H9_B[0-9]+_(D[0-9]+)$"
} else {
  "^IMR90_B[0-9]+_(D[0-9]+)$"
}

###############################################################
# Find complexes using EBI annotated list and prune complexes #
###############################################################

# Lookup table for EBI complexes
EBI_LT <- readRDS("00-LookupTables/ebi_cp_lookup.rds")

# Find complexes using EBI annotated list
complexes_EBI <- find_complexes(
  data_list = data_list,
  lookup_table = EBI_LT,
  id_col_data = "Accession",
  id_col_lookup = "uniprot_id"
)

# Build universe from all proteins identified across datasets (canonical accessions)
universe <- unique(sub(
  "-1$",
  "",
  trimws(unlist(lapply(data_list, function(dt) {
    unlist(lapply(dt[["Accession"]], function(x) {
      if (is.na(x) || x == "") {
        return(character(0))
      }
      trimws(strsplit(x, ";")[[1]])
    }))
  })))
))
universe <- universe[nzchar(universe) & !is.na(universe)]

# Prune complexes using statistical test
complexes_EBI <- prune_complexes_stat(
  all_complexes = complexes_EBI$complexes,
  lookup_table = EBI_LT,
  id_col_data = "Accession",
  id_col_lookup = "uniprot_id",
  universe = universe,
  min_found = 2,
  min_total = 2,
  alpha = 0.05,
  adjust_method = "BH"
)

# Cleanup
rm(universe)

###############################
# Remove NM data if requested #
###############################
if (!include_NM) {
  for (complex in names(complexes_EBI$keep)) {
    data <- complexes_EBI$keep[[complex]]
    if (!is.null(data)) {
      data <- data[datatype != "NM"]
      complexes_EBI$keep[[complex]] <- data
    }
  }
}

#####################################################
# Check for duplicate complexes in complexes vector #
#####################################################
if (length(complexes) != length(unique(complexes))) {
  stop(
    "Duplicate complex names found in 'complexes' vector. Please ensure all complex names are unique."
  )
}

##########################################################################
# Check that all specified complexes are present in the pruned complexes #
##########################################################################
missing_complexes <- setdiff(complexes, names(complexes_EBI$keep))
if (length(missing_complexes) > 0) {
  warning(
    paste(
      "The following complexes were specified but not found in the pruned complexes:",
      paste(missing_complexes, collapse = ", ")
    )
  )
}


##########################################################################
# Resolve ambiguous peptides: if a peptide maps to multiple proteins,    #
# keep it only when exactly one of those proteins is a complex member,   #
# updating the accession to that member. Remove otherwise.               #
##########################################################################
for (complex in names(complexes_EBI$keep)) {
  data <- complexes_EBI$keep[[complex]]

  data_clean <- resolve_ambiguous_peptides(
    data,
    accession_col = "Accession",
    position_col = "Position in master protein",
    lookup_table = EBI_LT,
    id_col_lookup = "uniprot_id"
  )

  complexes_EBI$keep[[complex]] <- data_clean
}

##########################################
# Analyze each complex for complexoforms #
##########################################

complex_obj <- list(
  dCF = list(),
  constraints = list(),
  vsclust = list(),
  quantification = list()
)

for (complex in complexes) {
  if (verbose) {
    cat(paste(rep("#", 80), collapse = ""), "\n")
    cat("\nAnalyzing complex:", complex, "\n")
  }
  result <- list()
  data <- complexes_EBI$keep[[complex]]
  if (is.null(data)) {
    warning(paste("Complex not found:", complex))
    next
  }

  # Create result directory for complex
  complex_dir <- paste0("10-Results/", gsub("/", "_", gsub(" ", "_", complex)))
  if (!dir.exists(complex_dir)) {
    dir.create(complex_dir)
  }

  # Find number of complexoform groups
  res_with_dcf <- discover_complexoforms(
    complex_data = data,
    intensity_cols = intensity_cols,
    group_col = "Gene name",
    peptide_col = "Peptide",
    alpha = 0.05,
    feature_adjust_method = "BH",
    min_total_non_na_frac = 0.40,
    deep_split = 2,
    min_conditions = 2,
    min_reps_per_condition = 2,
    minClusterSize = 2,
    cond_regex = cond_regex,
    canonical_label = "dCF0",
    singleton_label = "dCF-1",
    return_condition_tests = FALSE,
    verbose = verbose
  )

  result$dCF <- res_with_dcf
  complex_obj$dCF[[complex]] <- res_with_dcf

  # Run clustering using VSClust (with constraints if included)
  nclust <- find_nclust(res_with_dcf$peptide_results)

  # Create constraints for clustering if requested
  restriction_matrix <- NULL
  if (include_constraints) {
    if (verbose) {
      cat("\nBuilding constraints for complex:", complex, "\n")
    }

    # Create cannot-link constraints based on Pearson correlation
    cannotlink_matrix <- create_cannot_link_Pearson(
      data = res_with_dcf$peptide_results,
      complex_name = complex,
      id_cols = c("Gene name", "Peptide"),
      intensity_cols = intensity_cols,
      n_rep = n_rep,
      n_cond = n_cond,
      min_shared_cond = 10,
      cannot_link_th = 0.5,
      complex_info = EBI_LT,
      verbose = verbose
    )

    # Convert to VSClust restriction matrix via initial unconstrained clustering
    constraint_result <- vsclust_to_restrictions(
      data = res_with_dcf$peptide_results,
      cannotlink_matrix = cannotlink_matrix,
      n_clusters = nclust,
      n_rep = n_rep,
      n_cond = n_cond,
      id_cols = c("Gene name", "Peptide"),
      value_cols = intensity_cols,
      verbose = verbose
    )
    restriction_matrix <- constraint_result$restriction_matrix
    result$constraints <- list(
      cannotlink = cannotlink_matrix,
      restriction = restriction_matrix,
      initial_clustering = constraint_result$initial_clustering
    )
    complex_obj$constraints[[complex]] <- result$constraints
  }

  vsclust_res <- vsclust_on_complex(
    data = res_with_dcf$peptide_results,
    id_cols = c("Gene name", "Peptide"),
    intensity_cols = intensity_cols,
    n_rep = n_rep,
    n_cond = n_cond,
    n_clusters = nclust,
    grouped_replicates = TRUE,
    restriction_matrix = restriction_matrix,
    sds_amplification_factor = 1
  )

  vsclust_res$original_data <- data.table::copy(res_with_dcf$peptide_results)
  result$vsclust <- vsclust_res
  complex_obj$vsclust[[complex]] <- vsclust_res

  # Quantify complexoform abundance from VSClust memberships
  quantification_res <- quantify_complexoform_abundance_from_vsclust(
    vsclust_result = vsclust_res,
    intensity_cols = intensity_cols,
    id_cols = c("Gene name", "Peptide"),
    aggregation = "weighted_mean",
    membership_threshold = 0,
    normalize_membership = FALSE,
    cond_regex = cond_regex
  )
  result$quantification <- quantification_res
  complex_obj$quantification[[complex]] <- quantification_res

  # Save membership matrix summary
  cluster_summary <- make_cluster_summary(
    membership_matrix = vsclust_res$ClustOut$Bestcl$membership
  )
  data.table::fwrite(
    cluster_summary,
    file.path(complex_dir, "cluster_summary.csv")
  )

  # Save quantification outputs as tabular files for downstream analysis
  data.table::fwrite(
    quantification_res$abundance_wide,
    file.path(complex_dir, "complexoform_abundance.csv")
  )
  data.table::fwrite(
    quantification_res$group_stats,
    file.path(complex_dir, "complexoform_group_stats.csv")
  )

  # Save combined quantification plot (line + stacked bar, 16 x 6)
  save_complexoform_quantification_plot(
    quantification_result = quantification_res,
    file_path = file.path(complex_dir, "complexoform_quantification_plot.pdf"),
    complex_name = complex
  )

  # Save entire result object for complex
  saveRDS(
    result,
    file.path(complex_dir, "complexoform_results.rds")
  )

  # Create and save expression-membership plot for complex
  if (verbose) {
    cat("\nCreating visualization for complex:", complex, "\n")
  }
  grDevices::pdf(
    file = file.path(complex_dir, "complexoform_plot.pdf"),
    width = 16,
    height = 9,
    onefile = FALSE
  )
  create_complexoform_plot_from_pipeline(
    vsclust_result = vsclust_res,
    complex_name = complex,
    n_rep = n_rep,
    n_cond = n_cond
  )
  grDevices::dev.off()

  # Create and save datatype sequence plot for complex
  save_datatype_sequence_plot_pdf(
    complex_data = data,
    file_path = file.path(complex_dir, "datatype_sequence_plot.pdf"),
    complex_name = complex,
    ncol = 1
  )

  # Create and save per-cluster sequence plots (one PDF per cluster)
  save_cluster_sequence_plots_pdf(
    vsclust_result = vsclust_res,
    complex_data = data,
    output_dir = complex_dir,
    id_cols = c("Gene name", "Peptide"),
    ncol = 1
  )

  # Create and save PPCA plot colored by protein of origin and dCF grouping
  save_grouping_pca_pdf(
    vsclust_result = vsclust_res,
    intensity_cols = intensity_cols,
    n_rep = n_rep,
    cond_regex = cond_regex,
    complex_name = complex,
    file_path = file.path(complex_dir, "grouping_pca_plot.pdf")
  )
}
