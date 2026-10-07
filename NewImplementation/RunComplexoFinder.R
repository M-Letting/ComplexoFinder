# Set working directory to project root
setwd(here::here())

# Source functions
source("NewImplementation/00-LookupTables/GetLookupTables.R")
source("NewImplementation/02-Identification/IdentifyComplexes.R")
source("NewImplementation/02-Identification/DiscoverComplexoforms.R")
source("NewImplementation/04-Clustering/VsClust.R")
source("NewImplementation/04-Clustering/ClusterSummary.R")
source("NewImplementation/03-Constraints/Constraints.R")
source("NewImplementation/05-Quantification/Quantification.R")
source("NewImplementation/09-Visualization/ExpressionMembershipPlot.R")
source("NewImplementation/09-Visualization/SequencePlot.R")
source("NewImplementation/09-Visualization/QuantificationPlot.R")
source("NewImplementation/09-Visualization/GroupingProteinPCA.R")

################################################################################
## Options #####################################################################
################################################################################

# Choose which dataset to analyze
dataset <- "H9" # Options: "H9", "IMR90"

# Verbose output
verbose <- TRUE # TRUE for detailed output, FALSE for minimal output

# Include NM (protein level abundance) in analysis
include_NM <- FALSE # TRUE to include protein-level data, FALSE to exclude

# Specify complexes to analyze - for all options, see EBI complex lookup table
complexes <- c(
  "B-WICH chromatin remodelling complex",
  "CRD-mediated mRNA stability complex",
  "Multiaminoacyl-tRNA synthetase complex",
  "Major Spliceosomal B complex",
  "26S proteasome complex",
  "60S cytosolic large ribosomal subunit",
  "Nuclear pore complex",
  "Eukaryotic translation initiation factor 3 complex",
  "Dynein-1 complex, variant 1",
  "Intraflagellar transport complex B",
  "Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  "Neuron-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  "Neural progenitor-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
  "Calcineurin-Calmodulin-AKAP5 complex, gamma-R1 variant",
  "SNARE complex STX4-SNAP29-SEC22b",
  "Cortical microtubule stabilization complex, KANK1 variant",
  "LIFT actin modulation complex",
  "Dynein-1 complex, variant 4",
  "WASH complex, variant WASHC1/WASHC2C",
  "Cav1.2 voltage-gated calcium channel complex, CACNA2D1-CACNB3 variant",
  "Dynactin complex",
  "Laminin-213 complex"
)

# ---- Pruning: which complexes are analysed as an assembly --------------------
min_found <- 2 # Minimum number of detected member proteins
min_fraction <- 0.75 # Minimum fraction of a complex's members detected

# ---- Discovery of discordant peptides (all pruned complexes, one call) -------
alpha <- 0.05
min_total_non_na_frac <- 0 # Completeness filter: the only one in the pipeline
min_conditions <- 2
min_reps_per_condition <- 1
core_z <- 1.5
tau_aware_core <- FALSE # TRUE: second core pass with tau^2 in its distance
correction <- "two_stage" # "two_stage": Bonferroni within complex, then BH
#                           "global_bh": BH across all peptides only
deep_split <- 2 # Dynamic tree cut that groups discordant peptides into dCFs:
#                 0-4, higher gives more, smaller groups
reuse_discovery <- TRUE # Reuse the saved discovery if data and settings match

# ---- Clustering --------------------------------------------------------------
# Minimum number of conditions a peptide must be observed in to be clustered.
# Discovery tests every peptide with a profile; a profile on a handful of
# conditions is too short to place in a cluster.
min_conditions_clustering <- 10
scaling <- "center" # "center": each peptide centred on its own mean (keeps
#                      amplitude, as discovery does); "standardize": also
#                      scaled to unit variance (shape only)
sds_source <- "complex" # "complex": limma within the complex (as VSClust does)
#                         "global": the discovery noise model (all complexes,
#                                   one prior per datatype)
seed <- 1 # Seed for VSClust's random starts

# ---- Constraints -------------------------------------------------------------
include_constraints <- TRUE # TRUE to include constraints, FALSE to exclude
# Overlapping peptides whose profiles differ may not share a cluster.
cannot_link_grouping <- TRUE # TRUE: pattern groups per position (overlapping
#   peptides are grouped by profile and linked across groups only); FALSE:
#   every overlapping pair tested on its own (the calibrated rule)
cannot_link_tolerance <- NULL # Difference between two profiles (log2, RMS)
#   that does not count as a difference. NULL: tau from discovery; 0: any
#   difference beyond noise
cannot_link_th <- NULL # Adjusted p-value for a link. NULL: 0.05
min_shared_cond <- 10

# ---- Quantification and plots ------------------------------------------------
membership_cutoff <- 0.5 # Membership needed to count as a member of a group

################################################################################
## Run ComplexoFinder ##########################################################
################################################################################

#############
# Load data #
#############

if (!dataset %in% c("H9", "IMR90")) {
  stop("Invalid dataset choice. Please choose 'H9' or 'IMR90'.")
}

data_file <- file.path("NewImplementation/00-Data", paste0(dataset, "_data.rds"))

# Check if data files exist, if not run data processing script
if (!file.exists(data_file)) {
  source("NewImplementation/01-DataProcessing/ProcessData.R")
}
data_list <- readRDS(data_file)

intensity_cols <- grep(
  paste0("^", dataset, "_B[0-9]+_D[0-9]+$"),
  names(data_list[[1]]),
  value = TRUE
)
cond_regex <- paste0("^", dataset, "_B[0-9]+_(D[0-9]+)$")
id_cols <- c("Gene name", "Peptide")

# Lookup table for EBI complexes (downloaded and built on first use)
EBI_LT <- get_ebi_lookup()

# Protein sequences for the sequence plots
sequence_lookup <- get_sequence_lookup()

# Check for duplicate complexes in complexes vector
if (length(complexes) != length(unique(complexes))) {
  stop(
    "Duplicate complex names found in 'complexes' vector. Please ensure all complex names are unique."
  )
}

# Output directory
results_dir <- file.path("NewImplementation/10-Results", dataset)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################
# Find complexes using EBI annotated list and prune complexes #
###############################################################

complexes_EBI <- find_complexes(
  data_list = data_list,
  lookup_table = EBI_LT,
  data_id_column = "Accession",
  lookup_id_column = "uniprot_id"
)

# Remove NM data if requested
if (!include_NM) {
  complexes_EBI <- complexes_EBI[datatype != "NM"]
}

##########################################################################
# Resolve ambiguous peptides: if a peptide maps to multiple proteins,    #
# keep it only when exactly one member of the complex is among them,     #
# updating the accession to that member. Remove otherwise.               #
##########################################################################
complexes_EBI <- resolve_ambiguous_peptides(
  complexes_EBI,
  lookup_table = EBI_LT,
  accession_column = "Accession",
  position_column = "Position in master protein",
  gene_column = "Gene name",
  lookup_id_column = "uniprot_id"
)

# Keep complexes with enough of their members detected. This comes after the
# two steps above, so a member counts as detected only when it has a peptide
# that is actually analysed.
complexes_EBI <- prune_complexes(
  complex_data = complexes_EBI,
  lookup_table = EBI_LT,
  data_id_column = "Accession",
  lookup_id_column = "uniprot_id",
  min_proteins_found = min_found,
  min_fraction_found = min_fraction
)
data.table::fwrite(
  complexes_EBI$summary,
  file.path(results_dir, "complex_pruning_summary.csv")
)

# Check that all specified complexes are present in the pruned complexes
missing_complexes <- setdiff(complexes, unique(complexes_EBI$keep$complex_name))
if (length(missing_complexes) > 0) {
  warning(
    paste(
      "The following complexes were specified but not found in the pruned complexes:",
      paste(missing_complexes, collapse = ", ")
    )
  )
}

# Give every peptide an identifier that is unique within its complex
complex_data <- add_peptide_ids(complexes_EBI$keep)

#########################################################################
# Discover discordant peptides in all pruned complexes, in one call:    #
# variance priors, dispersion and FDR are estimated across complexes    #
#########################################################################

discovery_file <- file.path(results_dir, "discovery.rds")
discovery_settings <- list(
  data = c(nrow(complex_data), file.mtime(data_file)),
  lookup = nrow(EBI_LT),
  pruning = c(min_found, min_fraction, include_NM),
  discovery = c(
    alpha,
    min_total_non_na_frac,
    min_conditions,
    min_reps_per_condition,
    core_z,
    tau_aware_core
  ),
  labelling = c(correction, "dynamic", deep_split)
)

discovery <- NULL
if (reuse_discovery && file.exists(discovery_file)) {
  saved <- readRDS(discovery_file)
  if (identical(saved$settings, discovery_settings)) {
    discovery <- saved
    if (verbose) cat("Reusing saved discovery:", discovery_file, "\n")
  }
  rm(saved)
}

if (is.null(discovery)) {
  if (verbose) {
    cat(
      "Discovering discordant peptides in",
      data.table::uniqueN(complex_data$complex_id),
      "complexes...\n"
    )
  }
  discovery <- discover_complexoforms(
    data = complex_data,
    intensity_columns = intensity_cols,
    assembly_column = "complex_name",
    peptide_column = "Peptide",
    condition_regex = cond_regex,
    alpha = alpha,
    min_fraction_observed = min_total_non_na_frac,
    min_conditions = min_conditions,
    min_replicates_per_condition = min_reps_per_condition,
    datatype_column = "datatype",
    core_cut_height = core_z,
    variance_aware_core = tau_aware_core,
    correction = correction,
    deep_split = deep_split,
    min_cluster_size = 2,
    canonical_label = "dCF0",
    singleton_label = "dCF-1",
    verbose = verbose
  )
  discovery$settings <- discovery_settings
  saveRDS(discovery, discovery_file)
}

peptide_results <- discovery$peptide_results
peptide_results[, n_cond_observed := rowSums(is.finite(discovery$model$cond_mean))]

# Peptide-level calls for all complexes (without intensities)
data.table::fwrite(
  peptide_results[, !intensity_cols, with = FALSE],
  file.path(results_dir, "discovery_peptides.csv")
)
data.table::fwrite(
  peptide_results[,
    .(
      peptides = .N,
      tested = sum(!is.na(p_profile)),
      discordant = sum(discordant),
      n_dCF = data.table::uniqueN(dCF[!dCF %in% c("dCF0", "dCF-1")])
    ),
    by = .(complex_id, complex_name)
  ],
  file.path(results_dir, "discovery_summary.csv")
)

##########################################
# Analyze each complex for complexoforms #
##########################################

complex_obj <- list(
  dCF = list(),
  constraints = list(),
  vsclust = list(),
  quantification = list()
)
run_summary <- list()

for (complex in complexes) {
  if (verbose) {
    cat(paste(rep("#", 80), collapse = ""), "\n")
    cat("\nAnalyzing complex:", complex, "\n")
  }
  result <- list()

  # Peptides of this complex with a profile long enough to cluster
  rows <- which(
    peptide_results$complex_name == complex &
      discovery$model$tested &
      peptide_results$n_cond_observed >= min_conditions_clustering
  )
  data <- peptide_results[rows]
  if (nrow(data) < 3) {
    warning(paste("Complex not found or too few peptides:", complex))
    next
  }

  # Create result directory for complex
  complex_dir <- file.path(results_dir, gsub("/", "_", gsub(" ", "_", complex)))
  if (!dir.exists(complex_dir)) {
    dir.create(complex_dir)
  }

  result$dCF <- data
  complex_obj$dCF[[complex]] <- data

  # Number of clusters: the canonical group plus one per dCF
  nclust <- find_nclust(data)

  # Standard deviations for VSClust's fuzzifier
  sds <- if (sds_source == "global") sqrt(discovery$model$s2[rows]) else NULL

  vsclust_args <- list(
    data = data,
    id_columns = id_cols,
    intensity_columns = intensity_cols,
    condition_regex = cond_regex,
    n_clusters = nclust,
    standard_deviations = sds,
    scaling = scaling,
    seed = seed
  )

  # Create constraints for clustering if requested
  vsclust_res <- NULL
  n_cannot_link <- NA_integer_
  violated_initial <- NA_real_
  violated_final <- NA_real_
  converged_initial <- NA
  if (include_constraints && nclust > 1) {
    if (verbose) {
      cat("\nBuilding constraints for complex:", complex, "\n")
    }

    # Cannot-link constraints: overlapping peptides with differing profiles
    cannotlink_matrix <- create_cannot_link(
      data = data,
      id_columns = id_cols,
      intensity_columns = intensity_cols,
      condition_regex = cond_regex,
      grouping = cannot_link_grouping,
      tolerance = cannot_link_tolerance,
      link_threshold = cannot_link_th,
      min_shared_conditions = min_shared_cond,
      noise_model = list(
        squared_standard_error = discovery$model$se2[rows, , drop = FALSE],
        between_peptide_variance = discovery$tau2
      ),
      verbose = verbose
    )
    n_cannot_link <- sum(cannotlink_matrix) / 2

    # Initial unconstrained clustering, turned into a restriction matrix
    constraint_result <- do.call(
      vsclust_to_restrictions,
      c(
        vsclust_args,
        list(cannotlink_matrix = cannotlink_matrix, verbose = verbose)
      )
    )
    initial_clustering <- constraint_result$initial_clustering
    converged_initial <- initial_clustering$ClustOut$converged
    violated_initial <- cannot_link_violations(
      initial_clustering$ClustOut$Bestcl$cluster,
      cannotlink_matrix
    )[["violated"]]

    result$constraints <- list(
      cannotlink = cannotlink_matrix,
      restriction = constraint_result$restriction_matrix,
      initial_clustering = initial_clustering
    )
    complex_obj$constraints[[complex]] <- result$constraints

    if (n_cannot_link == 0) {
      # Nothing to constrain: the initial clustering is the result
      vsclust_res <- initial_clustering
    } else {
      # Constrained clustering, started from the initial clustering's centers
      # so that the restrictions refer to the same clusters
      vsclust_res <- do.call(
        vsclust_on_complex,
        c(
          vsclust_args,
          list(
            restriction_matrix = constraint_result$restriction_matrix,
            initial_centers = constraint_result$initial_centers
          )
        )
      )
    }
    violated_final <- cannot_link_violations(
      vsclust_res$ClustOut$Bestcl$cluster,
      cannotlink_matrix
    )[["violated"]]
    if (verbose) {
      cat(
        "Cannot-link pairs sharing a cluster:",
        round(100 * violated_initial),
        "% before,",
        round(100 * violated_final),
        "% after constraints\n"
      )
    }
  } else {
    vsclust_res <- do.call(vsclust_on_complex, vsclust_args)
  }

  result$vsclust <- vsclust_res
  complex_obj$vsclust[[complex]] <- vsclust_res

  # Quantify complexoform abundance from VSClust memberships
  quantification_res <- quantify_complexoform_abundance_from_vsclust(
    vsclust_result = vsclust_res,
    intensity_columns = intensity_cols,
    id_columns = id_cols,
    aggregation = "weighted_mean",
    membership_threshold = membership_cutoff,
    normalize_membership = FALSE,
    center_peptides = TRUE,
    condition_regex = cond_regex
  )
  result$quantification <- quantification_res
  complex_obj$quantification[[complex]] <- quantification_res

  run_summary[[complex]] <- data.table::data.table(
    complex = complex,
    peptides = nrow(data),
    discordant = sum(data$discordant),
    n_clusters = nclust,
    cannot_link_pairs = n_cannot_link,
    violated_initial = violated_initial,
    violated_final = violated_final,
    members = sum(
      matrixStats::rowMaxs(vsclust_res$ClustOut$Bestcl$membership) >=
        membership_cutoff
    ),
    # FALSE if VSClust stopped at its iteration limit: in the unconstrained
    # run the restrictions are built from, or in the final run
    converged_initial = converged_initial,
    converged_final = vsclust_res$ClustOut$converged
  )

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
    membership_cutoff = membership_cutoff
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
  unlink(list.files(
    complex_dir,
    pattern = "^cluster_sequence_plot_[0-9]+\\.pdf$",
    full.names = TRUE
  ))
  save_cluster_sequence_plots_pdf(
    vsclust_result = vsclust_res,
    complex_data = data,
    output_dir = complex_dir,
    id_cols = id_cols,
    ncol = 1,
    membership_cutoff = membership_cutoff
  )

  # Create and save PPCA plot colored by protein of origin and dCF grouping
  save_grouping_pca_pdf(
    vsclust_result = vsclust_res,
    intensity_cols = intensity_cols,
    cond_regex = cond_regex,
    complex_name = complex,
    file_path = file.path(complex_dir, "grouping_pca_plot.pdf")
  )
}

# One line per analysed complex
if (length(run_summary) > 0) {
  data.table::fwrite(
    data.table::rbindlist(run_summary),
    file.path(results_dir, "run_summary.csv")
  )
}
