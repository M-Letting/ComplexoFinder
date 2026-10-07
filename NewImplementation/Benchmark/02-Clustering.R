# ============================================================================
# 02-Clustering.R
# ============================================================================
# The steps after discovery - number of clusters, constraints, VSClust and
# quantification - old against new, end to end, on ProteoMaker, where the
# simulation says which peptides belong together and what their profile is.
#
# A protein is the assembly, and the 25 proteins with most peptides (in
# HighNoNA) are clustered in every dataset.
#
# Configurations:
#
#   old, unconstrained       old discovery -> number of clusters -> VSClust
#   old (pipeline)           + Pearson cannot-links via vsclust_to_restrictions()
#   new, unconstrained       new discovery (whole dataset) -> number of
#                            clusters -> VSClust
#   new (pipeline)           + Pearson cannot-links, constrained run started
#                            from the unconstrained run's centers. (The
#                            pipeline's default cannot-link rule is the one of
#                            "new, group links".)
#   new, calibrated links    create_cannot_link(grouping = FALSE, tolerance = 0)
#   new, ZRMSD links         create_cannot_link(method = "zrmsd")
#   new, CCC links           create_cannot_link(method = "ccc")
#   new, tolerance links (0.25)   grouping = FALSE, tolerance = 0.25: pairs
#                            whose profiles differ by more than 0.25 log2
#   new, group links         grouping = TRUE, tolerance = 0: pattern groups
#                            per position
#   new, group links (0.25)  the same with a tolerance of 0.25 log2
#   new, scaling = standardize   VSClust on standardized instead of centered
#                            profiles
#   new, sds = global        VSClust's standard deviations from the discovery
#                            noise model
#   new, tau-aware core      discovery with variance_aware_core = TRUE
#   new, global BH           discovery with correction = "global_bh"
#   new, true k (reference)  the new pipeline given the true number of groups
#
# The "links" configurations differ in the rule given to create_cannot_link();
# scaling, sds, tau-aware core and global BH are switches of
# NewImplementation/RunComplexoFinder.R.
#
# Every clustering is scored against the truth groups (score_clustering()) and
# quantified four ways (score_quantification()).
#
# The clusterings are saved in results/clustering_runs.rds and scored from
# there, so changing a metric does not need the pipelines to run again, and a
# new configuration only runs itself. Set `rerun_clustering <- TRUE` (or
# delete the file) to run everything again.
#
# Writes results/clustering_metrics.csv (one row per dataset, protein and
# configuration).
#
# Run from the project root (~40 min; a few minutes when the clusterings
# exist):
#   Rscript NewImplementation/Benchmark/02-Clustering.R
# ============================================================================

suppressMessages(library(data.table))
setwd(here::here())
source("NewImplementation/Benchmark/BenchmarkData.R")
source("NewImplementation/Benchmark/BenchmarkFunctions.R")

results_dir <- "NewImplementation/Benchmark/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

rerun_clustering <- FALSE
runs_file <- file.path(results_dir, "clustering_runs.rds")

new_pipeline <- load_new_implementation()

n_proteins <- 25
membership_cutoff <- 0.5
# Expected profiles closer than this (log2, root mean square) count as one
# true pattern in the tolerance version of the potential forms
pattern_tolerance <- 0.25
# ProteoMaker has 10 conditions; the pipeline's H9 settings (cluster peptides
# seen in >= 10 of 28 conditions, cannot-links need 10 shared conditions) are
# scaled to the same share of the conditions.
min_conditions_clustering <- 4
min_shared_conditions <- 3

# The proteins with most peptides in HighNoNA
top_proteins <- load_proteomaker_dataset("HighNoNA")$data[, .N, by = Accession][
  order(-N)][seq_len(n_proteins), Accession]

# (saved run, clustering within the run, configuration label)
configurations <- list(
  c("old", "none", "old, unconstrained"),
  c("old", "pearson", "old (pipeline)"),
  c("new", "none", "new, unconstrained"),
  c("new", "pearson", "new (pipeline)"),
  c("new", "calibrated", "new, calibrated links"),
  c("new_links", "zrmsd", "new, ZRMSD links"),
  c("new_links", "ccc", "new, CCC links"),
  c("new_links2", "tolerance 0.25", "new, tolerance links (0.25)"),
  c("new_links2", "groups", "new, group links"),
  c("new_links2", "groups 0.25", "new, group links (0.25)"),
  c("new_standardize", "pearson", "new, scaling = standardize"),
  c("new_sds", "pearson", "new, sds = global"),
  c("new_tau", "pearson", "new, tau-aware core"),
  c("new_bh", "pearson", "new, global BH"),
  c("new_true_k", "pearson", "new, true k (reference)")
)
run_names <- unique(vapply(configurations, `[`, character(1), 1))

# ---------------------------------------------------------------------------
# 1. Run the pipelines (only what is not saved yet)
# ---------------------------------------------------------------------------
# all_runs[[dataset]][[protein]][[run]] is the result of cluster_assembly_old()
# or cluster_assembly_new()
all_runs <- if (file.exists(runs_file) && !rerun_clustering) readRDS(runs_file) else list()

for (dataset_name in proteomaker_datasets) {
  missing_runs <- function(protein) {
    setdiff(run_names, names(all_runs[[dataset_name]][[protein]]))
  }
  dataset <- NULL

  for (protein in top_proteins) {
    runs_to_do <- missing_runs(protein)
    if (length(runs_to_do) == 0) next

    if (is.null(dataset)) {
      dataset <- load_proteomaker_dataset(dataset_name)
      cat(sprintf("\n[%s] %d peptides\n", dataset_name, nrow(dataset$data)))
      # New discovery: once per dataset and setting, as the pipeline does over
      # all complexes
      discoveries <- list(
        default = run_discovery_new(new_pipeline, dataset)$discovery,
        variance_aware_core = run_discovery_new(
          new_pipeline, dataset, variance_aware_core = TRUE
        )$discovery,
        global_bh = run_discovery_new(
          new_pipeline, dataset, correction = "global_bh"
        )$discovery
      )
      if ("old" %in% unlist(lapply(top_proteins, missing_runs))) {
        old_pipeline <- load_old_implementation()
      }
    }
    if (!protein %in% dataset$data$Accession) next

    assembly_rows <- which(dataset$data$Accession == protein)
    # Truth groups of at least two peptides: what the cluster number estimates
    true_n_clusters <- max(
      sum(table(dataset$truth$truth_group[assembly_rows]) >= 2),
      1
    )
    run_new <- function(discovery, link_rules = "pearson", ...) {
      cluster_assembly_new(
        new_pipeline, dataset, discovery, assembly_rows,
        min_conditions_clustering = min_conditions_clustering,
        min_shared_conditions = min_shared_conditions,
        link_rules = link_rules,
        ...
      )
    }
    run_functions <- list(
      old = function() tryCatch(
        cluster_assembly_old(
          old_pipeline, dataset, dataset$data[assembly_rows],
          min_shared_conditions = min_shared_conditions
        ),
        error = function(e) {
          warning(dataset_name, " ", protein, " old: ", conditionMessage(e))
          NULL
        }
      ),
      new = function() run_new(
        discoveries$default, link_rules = c("pearson", "calibrated")
      ),
      new_links = function() run_new(
        discoveries$default, link_rules = c("zrmsd", "ccc")
      ),
      new_links2 = function() run_new(discoveries$default, link_rules = list(
        "tolerance 0.25" = list(method = "tolerance", tolerance = 0.25),
        "groups" = list(method = "groups", tolerance = 0),
        "groups 0.25" = list(method = "groups", tolerance = 0.25)
      )),
      new_standardize = function() run_new(
        discoveries$default, scaling = "standardize"
      ),
      new_sds = function() run_new(
        discoveries$default, standard_deviation_source = "global"
      ),
      new_tau = function() run_new(discoveries$variance_aware_core),
      new_bh = function() run_new(discoveries$global_bh),
      new_true_k = function() run_new(
        discoveries$default, n_clusters = true_n_clusters
      )
    )
    for (run_name in runs_to_do) {
      # list(NULL) keeps a failed run in the list, so it is not tried again
      all_runs[[dataset_name]][[protein]][run_name] <- list(run_functions[[run_name]]())
    }
    clusters_used <- function(run_name) {
      run <- all_runs[[dataset_name]][[protein]][[run_name]]
      if (is.null(run)) "-" else run$k
    }
    cat(sprintf(
      "  %s: %d peptides, k true %d | old %s, new %s, global BH %s\n",
      protein, length(assembly_rows), true_n_clusters,
      clusters_used("old"), clusters_used("new"), clusters_used("new_bh")
    ))
  }
  # Saved after every dataset, so a stopped run keeps what it has
  if (!is.null(dataset)) saveRDS(all_runs, runs_file)
}

# ---------------------------------------------------------------------------
# 2. Score the saved clusterings
# ---------------------------------------------------------------------------
metrics <- list()

for (dataset_name in names(all_runs)) {
  dataset <- load_proteomaker_dataset(dataset_name)
  intensity_columns <- dataset$intensity_columns
  all_peptide_ids <- paste(
    dataset$data[[dataset$assembly_column]],
    dataset$data[[dataset$peptide_column]],
    sep = "_"
  )

  for (protein in names(all_runs[[dataset_name]])) {
    assembly_rows <- which(dataset$data$Accession == protein)
    assembly_data <- dataset$data[assembly_rows]
    residue_range <- tstrsplit(as.character(assembly_data$Position), "-", fixed = TRUE)
    truth <- data.table(
      id = all_peptide_ids[assembly_rows],
      truth_group = dataset$truth$truth_group[assembly_rows],
      pf1 = grepl("(^|\\|)1(\\||$)", assembly_data$Proteoform_ID),
      # For the potential forms
      start = suppressWarnings(as.numeric(residue_range[[1]])),
      stop = suppressWarnings(as.numeric(residue_range[[2]])),
      proteoforms = strsplit(as.character(assembly_data$Proteoform_ID), "|", fixed = TRUE)
    )
    expected_profiles <- dataset$expected[assembly_rows, , drop = FALSE]
    # True patterns up to a tolerance: peptides whose expected profiles are
    # all within `pattern_tolerance` log2 (root mean square over conditions)
    # of one another
    centered_expected <- expected_profiles - rowMeans(expected_profiles)
    truth[, pattern_tol := if (nrow(centered_expected) > 1) {
      stats::cutree(
        stats::hclust(
          stats::dist(centered_expected) / sqrt(ncol(centered_expected)),
          method = "complete"
        ),
        h = pattern_tolerance
      )
    } else {
      1L
    }]
    quantification_data <- data.table(
      Identifier = all_peptide_ids[assembly_rows],
      assembly_data[, ..intensity_columns]
    )
    # Truth groups of at least two peptides: what the cluster number estimates
    true_n_clusters <- max(sum(table(truth$truth_group) >= 2), 1)

    for (configuration in configurations) {
      run <- all_runs[[dataset_name]][[protein]][[configuration[1]]]
      if (is.null(run)) next
      membership <- run$clusterings[[configuration[2]]]
      metrics[[length(metrics) + 1]] <- cbind(
        data.table(
          dataset = dataset_name, protein = protein, config = configuration[3],
          peptides = length(assembly_rows), k = run$k, k_true = true_n_clusters
        ),
        score_clustering(
          membership, truth, expected_profiles, membership_cutoff,
          cannot_link = run$cannot_link[[configuration[2]]]
        ),
        score_quantification(
          new_pipeline, membership, quantification_data, intensity_columns,
          dataset$condition_regex, truth, expected_profiles, membership_cutoff
        )
      )
    }
  }
}

fwrite(rbindlist(metrics, fill = TRUE), file.path(results_dir, "clustering_metrics.csv"))
cat("\nWrote", file.path(results_dir, "clustering_metrics.csv"), "\n")
