# ============================================================================
# BENCHMARK FUNCTIONS
# ============================================================================
# Shared by the benchmark scripts in this folder:
#
#   - the two pipelines, each loaded into its own environment so that
#     functions with the same name (discover_complexoforms(),
#     vsclust_on_complex(), ...) do not overwrite each other:
#       old pipeline = ComplexoFinder/
#       new pipeline = NewImplementation/
#   - wrappers that run one step of each pipeline the way the pipeline runs it
#     and return the result in one common layout
#   - the metrics
#
# `dataset` is a list from a loader in BenchmarkData.R.
#
# Paths are relative to the project root.
# ============================================================================

suppressMessages({
  library(data.table)
  library(matrixStats)
})

# ---------------------------------------------------------------------------
# Pipelines
# ---------------------------------------------------------------------------

#' Source R files into an environment of their own
#'
#' @param files Character vector of file paths.
#' @return The environment holding everything the files define.
load_implementation <- function(files) {
  pipeline <- new.env(parent = globalenv())
  for (file in files) {
    suppressMessages(suppressWarnings(sys.source(file, envir = pipeline)))
  }
  pipeline
}

#' Load the old pipeline (ComplexoFinder/)
#'
#' @return Environment with the old pipeline's functions.
load_old_implementation <- function() {
  load_implementation(file.path("ComplexoFinder", c(
    "02-Identification/AssignComplexes.R",
    "02-Identification/IdentifyComplexes.R",
    "03-Constraints/Constraints.R",
    "04-Clustering/VsClustWrappers.R",
    "04-Clustering/VsClust.R",
    "05-Quantification/Quantification.R"
  )))
}

#' Load the new pipeline (NewImplementation/)
#'
#' @return Environment with the new pipeline's functions.
load_new_implementation <- function() {
  load_implementation(file.path("NewImplementation", c(
    "00-LookupTables/GetLookupTables.R",
    "02-Identification/IdentifyComplexes.R",
    "02-Identification/DiscoverComplexoforms.R",
    "04-Clustering/VsClust.R",
    "03-Constraints/Constraints.R",
    "05-Quantification/Quantification.R"
  )))
}

#' Run an expression without the cat(), message() and warning() output
#'
#' @param code An R expression.
#' @return The value of the expression.
quietly <- function(code) {
  value <- NULL
  invisible(utils::capture.output(
    value <- suppressMessages(suppressWarnings(code))
  ))
  value
}

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------
# Both wrappers return the calls in one layout, one row per input peptide:
#   group_id, peptide_id, score (the p-value to rank on; NA if not tested),
#   called (logical), dCF

#' Run the old pipeline's discovery on a dataset
#'
#' Calls the old `discover_complexoforms()` once per assembly, with the old
#' pipeline's settings.
#'
#' @param old_pipeline Environment from `load_old_implementation()`.
#' @param dataset Dataset list from a loader.
#'
#' @return data.table of calls: `group_id`, `peptide_id`, `score`, `called`,
#'   `dCF`. Assemblies with fewer than two peptides, or in which the old
#'   discovery fails, are left out.
run_discovery_old <- function(old_pipeline, dataset) {
  columns <- c(
    dataset$assembly_column,
    dataset$peptide_column,
    dataset$intensity_columns
  )
  assemblies <- split(
    dataset$data[, ..columns],
    dataset$data[[dataset$assembly_column]]
  )

  calls_by_assembly <- lapply(assemblies, function(assembly_data) {
    if (nrow(assembly_data) < 2) return(NULL)
    discovery <- tryCatch(
      quietly(old_pipeline$discover_complexoforms(
        complex_data = assembly_data,
        intensity_cols = dataset$intensity_columns,
        group_col = dataset$assembly_column,
        peptide_col = dataset$peptide_column,
        alpha = 0.05,
        feature_adjust_method = "BH",
        min_total_non_na_frac = 0.40,
        deep_split = 2,
        min_conditions = 2,
        min_reps_per_condition = 2,
        minClusterSize = 2,
        cond_regex = dataset$condition_regex,
        canonical_label = "dCF0",
        singleton_label = "dCF-1"
      )),
      error = function(e) NULL
    )
    if (is.null(discovery)) return(NULL)
    peptide_results <- discovery$peptide_results
    data.table(
      group_id = peptide_results[[dataset$assembly_column]],
      peptide_id = peptide_results[[dataset$peptide_column]],
      score = peptide_results$p_min_bh_feature,
      called = peptide_results$discordant,
      dCF = peptide_results$dCF
    )
  })
  rbindlist(calls_by_assembly)
}

#' Run the new pipeline's discovery on a dataset
#'
#' One `discover_complexoforms()` call over the whole dataset.
#'
#' @param new_pipeline Environment from `load_new_implementation()`.
#' @param dataset Dataset list from a loader.
#' @param ... Further arguments for `discover_complexoforms()`.
#'
#' @return List with `calls` (data.table: `group_id`, `peptide_id`, `score`,
#'   `called`, `dCF`) and `discovery` (the full result of
#'   `discover_complexoforms()`, with the noise model).
run_discovery_new <- function(new_pipeline, dataset, ...) {
  discovery <- new_pipeline$discover_complexoforms(
    data = dataset$data,
    intensity_columns = dataset$intensity_columns,
    assembly_column = dataset$assembly_column,
    peptide_column = dataset$peptide_column,
    condition_regex = dataset$condition_regex,
    datatype_column = dataset$datatype_column,
    ...
  )
  peptide_results <- discovery$peptide_results
  list(
    calls = data.table(
      group_id = peptide_results[[dataset$assembly_column]],
      peptide_id = peptide_results[[dataset$peptide_column]],
      score = peptide_results$p_profile,
      called = peptide_results$discordant,
      dCF = peptide_results$dCF
    ),
    discovery = discovery
  )
}

#' Area under the ROC curve, from ranks
#'
#' @param scores Numeric vector; a small score means "more likely positive".
#' @param is_positive Logical vector of the same length.
#' @return The AUROC, or `NA` if there are no positives or no negatives.
auc_rank <- function(scores, is_positive) {
  n_positive <- sum(is_positive)
  n_negative <- sum(!is_positive)
  if (n_positive == 0 || n_negative == 0) return(NA_real_)
  ranks <- rank(-scores)
  (sum(ranks[is_positive]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

#' Discovery metrics
#'
#' Every method is judged on the same peptides: every peptide of an assembly
#' with at least two peptides. A peptide a method did not score gets a p-value
#' of 1 and counts as not called.
#'
#' @param calls data.table from a discovery wrapper.
#' @param truth data.table with `group_id`, `peptide_id` and
#'   `truth_discordant`.
#' @return One-row data.table:
#'   peptides  number of peptides judged
#'   coverage  share of the peptides the method scored
#'   auroc     AUROC over all peptides
#'   within    AUROC inside each assembly holding both a discordant and a
#'             concordant peptide, pooled: how well the method tells a
#'             discordant peptide from the others of its assembly
#'   fpr       share called among peptides of assemblies with no discordant
#'             member (false alarms)
#'   tpr       share of discordant peptides called (power)
#'   fdp       share of the called peptides that are not discordant (false
#'             discovery proportion)
evaluate_discovery <- function(calls, truth) {
  peptides <- merge(truth, calls, by = c("group_id", "peptide_id"), all.x = TRUE)
  peptides[, assembly_size := .N, by = group_id]
  peptides <- peptides[assembly_size >= 2]
  peptides[, ranking_score := fifelse(is.finite(score), score, 1)]
  peptides[, called := !is.na(called) & called]
  peptides[, discordant_share := mean(truth_discordant), by = group_id]
  # AUROC within each mixed assembly, weighted by its number of
  # discordant-concordant pairs
  within_assembly <- peptides[discordant_share > 0 & discordant_share < 1, {
    n_discordant <- sum(truth_discordant)
    n_concordant <- sum(!truth_discordant)
    ranks <- rank(-ranking_score)
    list(
      auroc = (sum(ranks[truth_discordant]) - n_discordant * (n_discordant + 1) / 2) /
        (n_discordant * n_concordant),
      weight = as.numeric(n_discordant * n_concordant)
    )
  }, by = group_id]
  no_discordant <- peptides$discordant_share == 0
  data.table(
    peptides = nrow(peptides),
    coverage = mean(is.finite(peptides$score)),
    auroc = auc_rank(peptides$ranking_score, peptides$truth_discordant),
    within = if (nrow(within_assembly)) {
      sum(within_assembly$auroc * within_assembly$weight) / sum(within_assembly$weight)
    } else {
      NA_real_
    },
    fpr = if (any(no_discordant)) mean(peptides$called[no_discordant]) else NA_real_,
    tpr = if (any(peptides$truth_discordant)) {
      mean(peptides$called[peptides$truth_discordant])
    } else {
      NA_real_
    },
    fdp = if (any(peptides$called)) {
      mean(!peptides$truth_discordant[peptides$called])
    } else {
      NA_real_
    }
  )
}

# ---------------------------------------------------------------------------
# Clustering (one assembly at a time)
# ---------------------------------------------------------------------------
# Both wrappers return a list with
#   k            number of clusters VSClust was run with
#   clusterings  named list of membership matrices (rows = "<assembly>_<peptide>")
#   cannot_link  named list of the cannot-link matrices the constrained
#                clusterings were built from

#' Membership matrix that puts every peptide in one cluster
#'
#' @param peptide_ids Character vector of peptide identifiers.
#' @return One-column matrix of ones, with the identifiers as row names.
trivial_membership <- function(peptide_ids) {
  matrix(
    1, length(peptide_ids), 1,
    dimnames = list(peptide_ids, "membership of cluster 1")
  )
}

#' Run the old pipeline's clustering on one assembly
#'
#' Old discovery -> number of clusters -> Pearson cannot-links -> restrictions
#' from an unconstrained VSClust run -> constrained VSClust run, with the old
#' pipeline's settings.
#'
#' @param old_pipeline Environment from `load_old_implementation()`.
#' @param dataset Dataset list from a loader.
#' @param assembly_data data.table: the assembly's rows of `dataset$data`.
#' @param cannot_link_threshold Numeric. Pearson correlation below which two
#'   overlapping peptides are linked.
#' @param min_shared_conditions Integer. Conditions two peptides must share to
#'   be linked.
#'
#' @return List with `k`, `clusterings` (`none`: the unconstrained run,
#'   `pearson`: the constrained run) and `cannot_link` (`pearson`).
cluster_assembly_old <- function(
  old_pipeline,
  dataset,
  assembly_data,
  cannot_link_threshold = 0.5,
  min_shared_conditions = 3
) {
  intensity_columns <- dataset$intensity_columns
  column_conditions <- sub(dataset$condition_regex, "\\1", intensity_columns)
  n_conditions <- length(unique(column_conditions))
  n_replicates <- length(intensity_columns) / n_conditions
  peptide_ids <- paste(
    assembly_data[[dataset$assembly_column]],
    assembly_data[[dataset$peptide_column]],
    sep = "_"
  )

  discovery_results <- quietly(old_pipeline$discover_complexoforms(
    complex_data = assembly_data[,
      c(dataset$assembly_column, dataset$peptide_column, intensity_columns),
      with = FALSE
    ],
    intensity_cols = intensity_columns,
    group_col = dataset$assembly_column,
    peptide_col = dataset$peptide_column,
    alpha = 0.05, feature_adjust_method = "BH", min_total_non_na_frac = 0.40,
    deep_split = 2, min_conditions = 2, min_reps_per_condition = 2, minClusterSize = 2,
    cond_regex = dataset$condition_regex,
    canonical_label = "dCF0", singleton_label = "dCF-1"
  ))$peptide_results
  n_clusters <- old_pipeline$find_nclust(discovery_results)

  if (n_clusters == 1) {
    single_cluster <- trivial_membership(peptide_ids)
    return(list(
      k = 1L,
      clusterings = list(none = single_cluster, pearson = single_cluster),
      cannot_link = list()
    ))
  }

  # Cannot-link matrix, renamed to the identifiers the clustering uses
  link_id_columns <- c("Peptide", "Accession", "Position", "Proteoform_ID", "Peptidoform")
  link_data <- copy(assembly_data[, c(link_id_columns, intensity_columns), with = FALSE])
  link_ids <- do.call(paste, c(link_data[, ..link_id_columns], sep = "_"))
  cannot_link <- quietly(old_pipeline$create_cannot_link_Pearson(
    data = link_data, complex_name = NULL, id_cols = link_id_columns,
    identifier_info = link_id_columns, intensity_cols = intensity_columns,
    n_rep = n_replicates, n_cond = n_conditions,
    min_shared_cond = min_shared_conditions, cannot_link_th = cannot_link_threshold,
    complex_info = NULL
  ))
  dimnames(cannot_link) <- rep(
    list(peptide_ids[match(rownames(cannot_link), link_ids)]),
    2
  )

  restrictions <- quietly(old_pipeline$vsclust_to_restrictions(
    data = copy(discovery_results), cannotlink_matrix = cannot_link,
    n_clusters = n_clusters, n_rep = n_replicates, n_cond = n_conditions,
    id_cols = c(dataset$assembly_column, dataset$peptide_column),
    value_cols = intensity_columns
  ))

  clustering_data <- data.table(
    Identifier = peptide_ids,
    assembly_data[, ..intensity_columns]
  )
  constrained <- quietly(old_pipeline$vsclust_on_complex(
    data = clustering_data, id_cols = "Identifier", intensity_cols = intensity_columns,
    n_rep = n_replicates, n_cond = n_conditions, n_clusters = n_clusters,
    grouped_replicates = TRUE, restriction_matrix = restrictions$restriction_matrix
  ))

  list(
    k = n_clusters,
    clusterings = list(
      none = restrictions$initial_clustering$ClustOut$Bestcl$membership,
      pearson = constrained$ClustOut$Bestcl$membership
    ),
    cannot_link = list(pearson = cannot_link)
  )
}

#' Run the new pipeline's clustering on one assembly
#'
#' Given the dataset-wide discovery result: number of clusters from the dCF
#' labels -> unconstrained VSClust run -> cannot-links -> restrictions ->
#' constrained run started from the unconstrained run's centers. The
#' unconstrained run is shared by all cannot-link rules.
#'
#' @param new_pipeline Environment from `load_new_implementation()`.
#' @param dataset Dataset list from a loader.
#' @param discovery Result of the new `discover_complexoforms()` on the whole
#'   dataset.
#' @param assembly_rows Row numbers of the assembly in `dataset$data` (and in
#'   `discovery`).
#' @param min_conditions_clustering Peptides observed in fewer conditions are
#'   not clustered.
#' @param min_shared_conditions Integer. Conditions two peptides must share to
#'   be linked.
#' @param link_rules Cannot-link rules to build a constrained clustering for:
#'   method names, or a named list of argument lists for
#'   `create_cannot_link()` (e.g. `list("groups 0.25" = list(method =
#'   "groups", tolerance = 0.25))`); the clusterings are returned under these
#'   names.
#' @param scaling Character. Scaling of the profiles VSClust clusters.
#' @param standard_deviation_source Character. "complex": VSClust's standard
#'   deviations are estimated from the assembly; "global": they come from the
#'   discovery's noise model.
#' @param seed Integer. Seed of VSClust's random starts.
#' @param n_clusters Number of clusters to use instead of the one the dCF
#'   labels give.
#'
#' @return List with `k`, `clusterings` (`none`: the unconstrained run, and
#'   one per cannot-link rule) and `cannot_link` (one matrix per rule), or
#'   `NULL` when fewer than three peptides can be clustered.
cluster_assembly_new <- function(
  new_pipeline,
  dataset,
  discovery,
  assembly_rows,
  min_conditions_clustering,
  min_shared_conditions = 3,
  link_rules = c("pearson", "calibrated"),
  scaling = "center",
  standard_deviation_source = "complex",
  seed = 1,
  n_clusters = NULL
) {
  if (is.character(link_rules)) {
    link_rules <- setNames(
      lapply(link_rules, function(method) list(method = method)),
      link_rules
    )
  }
  intensity_columns <- dataset$intensity_columns
  id_columns <- c(dataset$assembly_column, dataset$peptide_column)
  conditions_observed <- rowSums(
    is.finite(discovery$model$cond_mean[assembly_rows, , drop = FALSE])
  )
  clustered_rows <- assembly_rows[
    discovery$model$tested[assembly_rows] &
      conditions_observed >= min_conditions_clustering
  ]
  clustered_data <- discovery$peptide_results[clustered_rows]
  peptide_ids <- do.call(paste, c(clustered_data[, ..id_columns], sep = "_"))
  if (nrow(clustered_data) < 3) {
    return(NULL)
  }

  if (is.null(n_clusters)) n_clusters <- new_pipeline$find_nclust(clustered_data)
  n_clusters <- min(n_clusters, nrow(clustered_data) - 1)
  if (n_clusters == 1) {
    single_cluster <- trivial_membership(peptide_ids)
    return(list(
      k = 1L,
      clusterings = c(
        list(none = single_cluster),
        setNames(rep(list(single_cluster), length(link_rules)), names(link_rules))
      ),
      cannot_link = list()
    ))
  }

  clustering_arguments <- list(
    data = clustered_data,
    id_columns = id_columns,
    intensity_columns = intensity_columns,
    condition_regex = dataset$condition_regex,
    n_clusters = n_clusters,
    scaling = scaling,
    seed = seed,
    cores = 2,
    standard_deviations = if (standard_deviation_source == "global") {
      sqrt(discovery$model$s2[clustered_rows])
    } else {
      NULL
    }
  )
  noise_model <- list(
    squared_standard_error = discovery$model$se2[clustered_rows, , drop = FALSE],
    between_peptide_variance = discovery$tau2
  )

  clusterings <- list()
  cannot_link_matrices <- list()
  initial_clustering <- NULL
  for (rule in names(link_rules)) {
    cannot_link <- do.call(new_pipeline$create_cannot_link, c(
      list(
        clustered_data,
        id_columns = id_columns,
        intensity_columns = intensity_columns,
        condition_regex = dataset$condition_regex,
        min_shared_conditions = min_shared_conditions,
        accession_column = dataset$assembly_column,
        position_column = "Position",
        noise_model = noise_model
      ),
      link_rules[[rule]]
    ))
    if (is.null(initial_clustering)) {
      # The unconstrained run is the same for every cannot-link rule
      restrictions <- do.call(
        new_pipeline$vsclust_to_restrictions,
        c(clustering_arguments, list(cannotlink_matrix = cannot_link))
      )
      initial_clustering <- restrictions$initial_clustering$ClustOut$Bestcl
      initial_centers <- restrictions$initial_centers
      clusterings$none <- initial_clustering$membership
      restriction_matrix <- restrictions$restriction_matrix
    } else {
      restriction_matrix <- restrictions_from_clustering(initial_clustering, cannot_link)
    }
    cannot_link_matrices[[rule]] <- cannot_link
    clusterings[[rule]] <- if (sum(cannot_link) == 0) {
      initial_clustering$membership
    } else {
      do.call(new_pipeline$vsclust_on_complex, c(
        clustering_arguments,
        list(restriction_matrix = restriction_matrix, initial_centers = initial_centers)
      ))$ClustOut$Bestcl$membership
    }
  }

  list(k = n_clusters, clusterings = clusterings, cannot_link = cannot_link_matrices)
}

#' Restriction matrix from a clustering and a cannot-link matrix
#'
#' Applies the rules of the new pipeline's `vsclust_to_restrictions()` to a
#' clustering that has already been run. For each cannot-link pair: peptides
#' in different clusters are each forbidden from the other's cluster; of two
#' peptides in one cluster, the one with the lower membership is forbidden
#' from it. A peptide forbidden from every cluster is allowed into all.
#'
#' @param best_clustering VSClust result (`ClustOut$Bestcl`) with `cluster`
#'   and `membership`.
#' @param cannot_link Peptide x peptide cannot-link matrix.
#'
#' @return Logical matrix, peptides x clusters; `TRUE` forbids the peptide
#'   from the cluster.
restrictions_from_clustering <- function(best_clustering, cannot_link) {
  cluster <- best_clustering$cluster
  n_clusters <- ncol(best_clustering$membership)
  restriction_matrix <- matrix(
    FALSE, length(cluster), n_clusters,
    dimnames = list(names(cluster), seq_len(n_clusters))
  )
  linked_pairs <- which(cannot_link & upper.tri(cannot_link), arr.ind = TRUE)
  for (pair in seq_len(nrow(linked_pairs))) {
    first_peptide <- rownames(cannot_link)[linked_pairs[pair, 1]]
    second_peptide <- rownames(cannot_link)[linked_pairs[pair, 2]]
    if (!first_peptide %in% names(cluster) || !second_peptide %in% names(cluster)) next
    first_cluster <- cluster[[first_peptide]]
    second_cluster <- cluster[[second_peptide]]
    if (first_cluster != second_cluster) {
      restriction_matrix[first_peptide, second_cluster] <- TRUE
      restriction_matrix[second_peptide, first_cluster] <- TRUE
    } else if (best_clustering$membership[first_peptide, first_cluster] >
      best_clustering$membership[second_peptide, first_cluster]) {
      restriction_matrix[second_peptide, first_cluster] <- TRUE
    } else {
      restriction_matrix[first_peptide, first_cluster] <- TRUE
    }
  }
  restriction_matrix[rowSums(restriction_matrix) == n_clusters, ] <- FALSE
  restriction_matrix
}

# ---------------------------------------------------------------------------
# Clustering metrics (ground truth: ProteoMaker)
# ---------------------------------------------------------------------------

#' Potential forms of a clustering
#'
#' A cluster is read as one form. Its member peptides are laid out along the
#' sequence and merged into stretches covered by peptides that overlap one
#' another. Every stretch can be filled in several ways, so the cluster stands
#' for
#'
#'   product over its stretches of (number of alternatives in the stretch)
#'
#' potential forms. The alternatives of a stretch are its member peptides
#' (`alternative_labels = NULL`), or the distinct labels those peptides carry
#' (e.g. their proteoform IDs). The clustering's value is the sum over
#' clusters, and within a cluster over proteins when `accession` is given.
#' Lower is better. A peptide that is no cluster's member does not count.
#'
#' With proteoform IDs as labels this is `PM_evaluate_total()` of
#' Benchmark/03-Clustering/05-CompareClustering.R, without labels
#' `PM_evaluate_position()`, and with `accession` `PC_evaluate_position()`.
#'
#' @param member Logical matrix, peptides x clusters: the peptide is a member
#'   of the cluster.
#' @param start,stop Residue range of each peptide (`NA`: the peptide is left
#'   out).
#' @param alternative_labels Optional list, one vector of labels per peptide.
#' @param accession Optional protein of each peptide.
#'
#' @return The number of potential forms.
potential_forms <- function(
  member,
  start,
  stop,
  alternative_labels = NULL,
  accession = NULL
) {
  if (is.null(accession)) accession <- rep("", nrow(member))
  has_position <- !is.na(start) & !is.na(stop) & !is.na(accession)
  total_forms <- 0
  for (cluster in seq_len(ncol(member))) {
    cluster_members <- which(member[, cluster] & has_position)
    for (protein in unique(accession[cluster_members])) {
      protein_members <- cluster_members[accession[cluster_members] == protein]
      protein_members <- protein_members[
        order(start[protein_members], stop[protein_members])
      ]
      # A new stretch starts where a peptide begins after everything before it
      # has ended
      stretch <- cumsum(c(
        TRUE,
        start[protein_members][-1] >
          cummax(stop[protein_members])[-length(protein_members)]
      ))
      alternatives <- if (is.null(alternative_labels)) {
        tabulate(stretch)
      } else {
        vapply(
          split(alternative_labels[protein_members], stretch),
          function(labels) length(unique(unlist(labels))),
          numeric(1)
        )
      }
      total_forms <- total_forms + prod(alternatives)
    }
  }
  total_forms
}

#' Adjusted Rand index between two labelings of the same items
#'
#' @param first_labels,second_labels Vectors of group labels, one per item.
#' @return The adjusted Rand index; 1 when both labelings put everything in
#'   one group or everything apart.
ari <- function(first_labels, second_labels) {
  contingency <- table(first_labels, second_labels)
  n_items <- sum(contingency)
  pairs_together_in_both <- sum(choose(contingency, 2))
  pairs_together_in_first <- sum(choose(rowSums(contingency), 2))
  pairs_together_in_second <- sum(choose(colSums(contingency), 2))
  expected_pairs <- pairs_together_in_first * pairs_together_in_second /
    choose(n_items, 2)
  maximum_above_expected <- (pairs_together_in_first + pairs_together_in_second) / 2 -
    expected_pairs
  if (maximum_above_expected == 0) {
    1
  } else {
    (pairs_together_in_both - expected_pairs) / maximum_above_expected
  }
}

#' Score one clustering of one assembly against the ground truth
#'
#' Every method is judged on all peptides of the assembly. A peptide is
#' assigned to the cluster of its highest membership if that membership is at
#' least `membership_cutoff`; peptides a method did not cluster, or that are a
#' member of no cluster, are unassigned.
#'
#' @param membership Membership matrix (rows named by peptide id).
#' @param truth data.table of the assembly's peptides: `id`, `truth_group`,
#'   `pf1` (peptide carries proteoform 1). For the potential forms also
#'   `start`, `stop` (residue range), `proteoforms` (list of proteoform IDs)
#'   and optionally `pattern_tol` (an id shared by peptides whose expected
#'   profiles are close enough to count as one pattern).
#' @param expected_profiles Numeric matrix, row-aligned with `truth`.
#' @param membership_cutoff Numeric. Membership needed to be a member of a
#'   cluster.
#' @param cannot_link Optional cannot-link matrix the clustering was built
#'   from.
#' @param min_form_size Integer. Peptides a cluster or a truth group needs to
#'   count as a form.
#'
#' @return One-row data.table:
#'   truth_groups  number of truth groups among the assembly's peptides
#'   clustered     share of the assembly's peptides the method clustered
#'   assigned      share of the assembly's peptides assigned to a cluster
#'   pair_precision  of the peptide pairs put in one cluster, the share that
#'                 belong to one truth group (low = clusters mix groups)
#'   pair_recall   of the peptide pairs in one truth group, the share put in
#'                 one cluster (low = groups are split or left unassigned)
#'   ari           adjusted Rand index against the truth groups, over all
#'                 peptides; each unassigned peptide is a group of its own.
#'                 With a single truth group it is 1 for one cluster holding
#'                 every peptide and 0 for anything else.
#'   ari_assigned  the same over the assigned peptides only
#'   purity        share of a cluster's assigned peptides that come from its
#'                 most common truth group, averaged over clusters weighted by
#'                 size
#'   profile_r2    share of the variation between the peptides' expected
#'                 profiles that the clusters explain (unassigned peptides
#'                 form one extra group); NA when all expected profiles are
#'                 equal
#'   top_fraction  of the peptides that carry proteoform 1, the largest share
#'                 assigned to one cluster
#'   potential_proteoforms, potential_position
#'                 potential forms (see potential_forms()) with the distinct
#'                 proteoform IDs, or the member peptides, of a stretch as its
#'                 alternatives. Lower is better.
#'   sum_unique_proteoforms
#'                 number of distinct proteoform IDs among a cluster's
#'                 members, summed over clusters. Lower is better.
#'   potential_proteoforms_1k, potential_position_1k, sum_unique_proteoforms_1k
#'                 the same three with a flexible cutoff: a peptide is a member
#'                 of every cluster it has a membership of at least 1/k in
#'                 (k = number of clusters). Every clustered peptide is then a
#'                 member of at least one cluster, and can be a member of
#'                 several.
#'   potential_proteoforms_top, potential_position_top
#'                 the same with every clustered peptide counted once, in the
#'                 cluster it has the highest membership in
#'   potential_pattern, potential_pattern_1k, potential_pattern_top
#'                 potential forms with the distinct truth groups (expected
#'                 profiles) of a stretch's member peptides as its
#'                 alternatives, for the three definitions of membership
#'   potential_pattern_tol_1k, potential_pattern_tol_top
#'                 the same with the `pattern_tol` ids of `truth` as the
#'                 alternatives
#'   forms_found   number of clusters with at least `min_form_size` members
#'   forms_recovered  number of truth groups of at least `min_form_size`
#'                 peptides that one cluster recovers: it holds more than half
#'                 of the group's peptides and more than half of its members
#'                 come from the group
#'   cl_pairs, cl_violated  cannot-link pairs, and the share of them that sit
#'                 in the same cluster
#'   The potential-form columns are NA unless `truth` has `start`, `stop` and
#'   `proteoforms`.
score_clustering <- function(
  membership,
  truth,
  expected_profiles,
  membership_cutoff = 0.5,
  cannot_link = NULL,
  min_form_size = 2
) {
  membership_row <- match(truth$id, rownames(membership))
  memberships <- membership[membership_row, , drop = FALSE]
  memberships[is.na(memberships)] <- 0
  top_cluster <- max.col(memberships, ties.method = "first")
  assigned <- rowMaxs(memberships) >= membership_cutoff
  cluster_label <- ifelse(
    assigned,
    as.character(top_cluster),
    paste0("unassigned_", seq_along(top_cluster))
  )

  # Pairs of peptides: together in a cluster / in a truth group / in both
  contingency <- table(cluster_label, truth$truth_group)
  pairs_in_both <- sum(choose(contingency, 2))
  pairs_in_cluster <- sum(choose(rowSums(contingency), 2))
  pairs_in_truth_group <- sum(choose(colSums(contingency), 2))

  purity <- if (any(assigned)) {
    assigned_contingency <- table(top_cluster[assigned], truth$truth_group[assigned])
    sum(apply(assigned_contingency, 1, max)) / sum(assigned_contingency)
  } else {
    NA_real_
  }

  # Variation between expected profiles (each centered on its own mean):
  # total, and what remains within the clusters
  centered_expected <- expected_profiles - rowMeans(expected_profiles)
  cluster_or_unassigned <- ifelse(assigned, as.character(top_cluster), "unassigned")
  sum_of_squares <- function(profiles) sum(sweep(profiles, 2, colMeans(profiles))^2)
  total_sum_of_squares <- sum_of_squares(centered_expected)
  within_sum_of_squares <- sum(vapply(
    unique(cluster_or_unassigned),
    function(group) {
      sum_of_squares(centered_expected[cluster_or_unassigned == group, , drop = FALSE])
    },
    numeric(1)
  ))

  top_fraction <- if (any(truth$pf1)) {
    max(colMeans(memberships[truth$pf1, , drop = FALSE] >= membership_cutoff))
  } else {
    NA_real_
  }

  # Forms: reported (clusters with members) and true (truth groups)
  truth_group_sizes <- table(truth$truth_group)
  cluster_sizes <- table(top_cluster[assigned])
  true_forms <- names(truth_group_sizes)[truth_group_sizes >= min_form_size]
  found_forms <- names(cluster_sizes)[cluster_sizes >= min_form_size]
  forms_recovered <- 0L
  if (length(true_forms) > 0 && length(found_forms) > 0) {
    overlap <- table(
      factor(top_cluster[assigned], levels = names(cluster_sizes)),
      factor(truth$truth_group[assigned], levels = names(truth_group_sizes))
    )[found_forms, true_forms, drop = FALSE]
    cluster_mostly_from_form <- overlap / as.vector(cluster_sizes[found_forms]) > 0.5
    cluster_holds_most_of_form <-
      sweep(overlap, 2, as.vector(truth_group_sizes[true_forms]), "/") > 0.5
    forms_recovered <- sum(
      colSums(cluster_mostly_from_form & cluster_holds_most_of_form) > 0
    )
  }

  # Potential forms, for three definitions of membership
  potential <- c(
    proteoforms = NA_real_, position = NA_real_, sum_unique = NA_real_,
    proteoforms_1k = NA_real_, position_1k = NA_real_, sum_unique_1k = NA_real_,
    proteoforms_top = NA_real_, position_top = NA_real_,
    pattern = NA_real_, pattern_1k = NA_real_, pattern_top = NA_real_,
    pattern_tol_1k = NA_real_, pattern_tol_top = NA_real_
  )
  if (all(c("start", "stop", "proteoforms") %in% names(truth))) {
    member <- memberships >= membership_cutoff
    member_flexible <- memberships >= 1 / ncol(memberships)
    member_top <- matrix(FALSE, nrow(memberships), ncol(memberships))
    clustered <- which(!is.na(membership_row))
    member_top[cbind(clustered, top_cluster[clustered])] <- TRUE
    count_unique_proteoforms <- function(member_matrix) {
      sum(apply(member_matrix, 2, function(is_member) {
        length(unique(unlist(truth$proteoforms[is_member])))
      }))
    }
    count_forms <- function(member_matrix, alternative_labels = NULL) {
      potential_forms(member_matrix, truth$start, truth$stop, alternative_labels)
    }
    truth_group_labels <- as.list(truth$truth_group)
    potential <- c(
      proteoforms = count_forms(member, truth$proteoforms),
      position = count_forms(member),
      sum_unique = count_unique_proteoforms(member),
      proteoforms_1k = count_forms(member_flexible, truth$proteoforms),
      position_1k = count_forms(member_flexible),
      sum_unique_1k = count_unique_proteoforms(member_flexible),
      proteoforms_top = count_forms(member_top, truth$proteoforms),
      position_top = count_forms(member_top),
      pattern = count_forms(member, truth_group_labels),
      pattern_1k = count_forms(member_flexible, truth_group_labels),
      pattern_top = count_forms(member_top, truth_group_labels),
      pattern_tol_1k = NA_real_,
      pattern_tol_top = NA_real_
    )
    if ("pattern_tol" %in% names(truth)) {
      pattern_labels <- as.list(truth$pattern_tol)
      potential[["pattern_tol_1k"]] <- count_forms(member_flexible, pattern_labels)
      potential[["pattern_tol_top"]] <- count_forms(member_top, pattern_labels)
    }
  }

  link_violations <- c(pairs = NA_real_, violated = NA_real_)
  if (!is.null(cannot_link)) {
    hard_cluster <- setNames(top_cluster, truth$id)
    hard_cluster[is.na(membership_row)] <- NA
    linked_pairs <- which(cannot_link & upper.tri(cannot_link), arr.ind = TRUE)
    first_cluster <- hard_cluster[rownames(cannot_link)[linked_pairs[, 1]]]
    second_cluster <- hard_cluster[rownames(cannot_link)[linked_pairs[, 2]]]
    both_clustered <- !is.na(first_cluster) & !is.na(second_cluster)
    link_violations <- c(
      pairs = sum(both_clustered),
      violated = if (any(both_clustered)) {
        mean(first_cluster[both_clustered] == second_cluster[both_clustered])
      } else {
        NA_real_
      }
    )
  }

  data.table(
    truth_groups = uniqueN(truth$truth_group),
    clustered = mean(!is.na(membership_row)),
    assigned = mean(assigned),
    pair_precision = if (pairs_in_cluster > 0) pairs_in_both / pairs_in_cluster else NA_real_,
    pair_recall = if (pairs_in_truth_group > 0) {
      pairs_in_both / pairs_in_truth_group
    } else {
      NA_real_
    },
    ari = ari(cluster_label, truth$truth_group),
    ari_assigned = if (sum(assigned) > 1) {
      ari(top_cluster[assigned], truth$truth_group[assigned])
    } else {
      NA_real_
    },
    purity = purity,
    profile_r2 = if (total_sum_of_squares > 1e-12) {
      1 - within_sum_of_squares / total_sum_of_squares
    } else {
      NA_real_
    },
    top_fraction = top_fraction,
    potential_proteoforms = potential[["proteoforms"]],
    potential_position = potential[["position"]],
    sum_unique_proteoforms = potential[["sum_unique"]],
    potential_proteoforms_1k = potential[["proteoforms_1k"]],
    potential_position_1k = potential[["position_1k"]],
    sum_unique_proteoforms_1k = potential[["sum_unique_1k"]],
    potential_proteoforms_top = potential[["proteoforms_top"]],
    potential_position_top = potential[["position_top"]],
    potential_pattern = potential[["pattern"]],
    potential_pattern_1k = potential[["pattern_1k"]],
    potential_pattern_top = potential[["pattern_top"]],
    potential_pattern_tol_1k = potential[["pattern_tol_1k"]],
    potential_pattern_tol_top = potential[["pattern_tol_top"]],
    forms_found = length(found_forms),
    forms_recovered = forms_recovered,
    cl_pairs = link_violations[["pairs"]],
    cl_violated = link_violations[["violated"]]
  )
}

#' Error of the quantified cluster profiles against the truth
#'
#' The simulation fixes how a group of peptides deviates from the rest of its
#' protein, not the profile itself: every protein also has a profile shared by
#' all its peptides, which is not in the metadata. So for every cluster with
#' members (membership >= `membership_cutoff`) the comparison is
#'
#'   quantified:  cluster profile - mean profile of all the protein's peptides
#'   truth:       mean expected profile of the members - mean expected profile
#'                of all the protein's peptides
#'
#' both centered across conditions. The result is the root mean squared
#' difference (log2), averaged over clusters.
#'
#' The same clustering is quantified four ways:
#'   who contributes   all:     every peptide, weighted by its membership
#'                              (the old pipeline's setting)
#'                     members: only peptides with membership >=
#'                              `membership_cutoff` (the new pipeline's)
#'   how               plain:   weighted mean of the intensities (old)
#'                     centred: each peptide centered on its own level first
#'                              (new)
#'
#' @param new_pipeline Environment from `load_new_implementation()`.
#' @param membership Membership matrix (rows named by peptide id).
#' @param data data.table with `Identifier` and the intensity columns: all
#'   peptides of the assembly, row-aligned with `truth` and
#'   `expected_profiles`.
#' @param intensity_columns Character vector of intensity column names.
#' @param condition_regex Regex with one capture group for the condition
#'   label.
#' @param truth data.table of the assembly's peptides, with `id`.
#' @param expected_profiles Numeric matrix, row-aligned with `truth`.
#' @param membership_cutoff Numeric. Membership needed to be a member of a
#'   cluster.
#'
#' @return One-row data.table: quant_all_plain, quant_all_centred,
#'   quant_members_plain, quant_members_centred.
score_quantification <- function(
  new_pipeline,
  membership,
  data,
  intensity_columns,
  condition_regex,
  truth,
  expected_profiles,
  membership_cutoff = 0.5
) {
  center_rows <- function(profiles) profiles - rowMeans(profiles, na.rm = TRUE)

  # Truth: the members' expected deviation from the protein
  truth_row <- match(rownames(membership), truth$id)
  member <- membership >= membership_cutoff
  true_deviation <- t(vapply(
    seq_len(ncol(membership)),
    function(cluster) {
      if (!any(member[, cluster])) return(rep(NA_real_, ncol(expected_profiles)))
      colMeans(expected_profiles[truth_row[member[, cluster]], , drop = FALSE]) -
        colMeans(expected_profiles)
    },
    numeric(ncol(expected_profiles))
  ))
  true_deviation <- center_rows(true_deviation)

  # The protein's mean profile, from all its peptides (each centered on its
  # own level)
  condition_mean <- new_pipeline$condition_means(
    as.matrix(data[, ..intensity_columns]),
    new_pipeline$parse_conditions(intensity_columns, condition_regex)
  )
  protein_profile <- colMeans(center_rows(condition_mean), na.rm = TRUE)

  quantification_error <- function(membership_threshold, center_peptides) {
    abundance <- new_pipeline$quantify_complexoform_abundance(
      peptide_data = data,
      membership_matrix = membership,
      intensity_columns = intensity_columns,
      aggregation = "weighted_mean",
      membership_threshold = membership_threshold,
      center_peptides = center_peptides,
      condition_regex = condition_regex
    )$abundance_wide
    quantified_deviation <- center_rows(
      sweep(center_rows(as.matrix(abundance[, -1])), 2, protein_profile)
    )
    mean(
      sqrt(rowMeans((quantified_deviation - true_deviation)^2, na.rm = TRUE)),
      na.rm = TRUE
    )
  }
  data.table(
    quant_all_plain = quantification_error(0, FALSE),
    quant_all_centred = quantification_error(0, TRUE),
    quant_members_plain = quantification_error(membership_cutoff, FALSE),
    quant_members_centred = quantification_error(membership_cutoff, TRUE)
  )
}
