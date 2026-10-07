# ============================================================================
# 03-RealData.R
# ============================================================================
# The two pipelines on the real H9 complexes. There is no ground truth here,
# so this is a description of what each pipeline does, plus the one thing that
# can be checked objectively: whether the cannot-link constraints are enforced.
#
#   old  ComplexoFinder/ run on its own data and lookup
#        (ComplexoFinder/00-Data/H9_data.rds, 00-LookupTables/ebi_cp_lookup.rds),
#        step by step as its RunComplexoFinder.R does
#   new  the saved results of NewImplementation/RunComplexoFinder.R
#        (NewImplementation/10-Results/H9/) - run that script first
#
# The complexes are those in the new pipeline's run summary.
#
# Per complex and pipeline (results/realdata_metrics.csv):
#   peptides          peptides that are clustered
#   discordant        share of them called discordant
#   k                 number of clusters
#   cl_pairs          cannot-link pairs
#   violated_initial  share of those pairs in one cluster, unconstrained run
#   violated_final    the same after the constrained run
#   members           share of peptides with a membership >= 0.5
#   nonmember_weight  share of a group's quantification weight that comes from
#                     peptides that are not its members (median over groups)
#
# Per complex, pipeline and restriction (results/realdata_restrictions.csv):
#   restriction  none (the unconstrained run), pearson (the old pipeline's
#                rule), and for the new pipeline: groups (its default rule:
#                pattern groups per position, default tolerance), tolerance
#                (pair by pair, default tolerance), calibrated (pair by pair,
#                no tolerance), groups without a tolerance, and the rules
#                pearson, zrmsd and ccc. The new pipeline's clusterings are
#                built from the saved unconstrained run and clustered as the
#                pipeline does it.
#   k            clusters asked for
#   cl_pairs     cannot-link pairs the rule gives
#   violated     share of them in one cluster
#   forms        clusters with at least two members: the complexoforms reported
#   members      share of peptides with a membership >= 0.5
#   potential_position      potential forms (potential_forms() in
#                BenchmarkFunctions.R) with the member peptides of a stretch as
#                its alternatives, counted per protein within each cluster.
#                Lower is better.
#   potential_position_1k   the same with a peptide being a member of every
#                cluster it has a membership >= 1/k in
#   potential_position_top  the same with every peptide counted in the cluster
#                it has the highest membership in
#
# Run from the project root (~15 min):
#   Rscript NewImplementation/Benchmark/03-RealData.R
# ============================================================================

suppressMessages(library(data.table))
setwd(here::here())
source("NewImplementation/Benchmark/BenchmarkFunctions.R")

results_dir <- "NewImplementation/Benchmark/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

new_results_dir <- "NewImplementation/10-Results/H9"
if (!file.exists(file.path(new_results_dir, "run_summary.csv"))) {
  stop("Run NewImplementation/RunComplexoFinder.R first (dataset H9).")
}
new_run_summary <- fread(file.path(new_results_dir, "run_summary.csv"))
complexes <- new_run_summary$complex
membership_cutoff <- 0.5

#' Share of cannot-link pairs in one cluster
#'
#' @param cluster Named vector of hard cluster assignments.
#' @param cannot_link Peptide x peptide cannot-link matrix.
#' @return The share, or `NA` when no linked pair has both peptides clustered.
violations <- function(cluster, cannot_link) {
  linked_pairs <- which(cannot_link & upper.tri(cannot_link), arr.ind = TRUE)
  first_cluster <- cluster[rownames(cannot_link)[linked_pairs[, 1]]]
  second_cluster <- cluster[rownames(cannot_link)[linked_pairs[, 2]]]
  both_clustered <- !is.na(first_cluster) & !is.na(second_cluster)
  if (any(both_clustered)) {
    mean(first_cluster[both_clustered] == second_cluster[both_clustered])
  } else {
    NA_real_
  }
}

#' Protein and residue range of clustered peptides, read from their identifiers
#'
#' The identifiers have the form "GENE_ACCESSION_ACCESSION [start-end]_...".
#' The range is read from the square brackets (an isoform accession such as
#' "P28066-2" also looks like a range).
#'
#' @param peptide_ids Character vector of peptide identifiers.
#'
#' @return data.table with `id`, `accession`, `start` and `stop`; `NA` for an
#'   identifier without an accession or without a range ("[NOT FOUND]").
peptide_positions <- function(peptide_ids) {
  # The identifier without the gene name
  label <- sub("^[^_]+_", "", peptide_ids)
  accession_match <- regmatches(label, regexpr(
    "([OPQ][0-9][A-Z0-9]{3}[0-9](?:-[0-9]+)?|[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-[0-9]+)?)",
    label, perl = TRUE
  ), invert = NA)
  # Regex match of each label: the whole match, the start and the end
  range_match <- regmatches(label, regexec("\\[(\\d+)\\s*[-\u2013]\\s*(\\d+)\\]", label))
  matched_number <- function(match_index) {
    suppressWarnings(as.numeric(vapply(
      range_match,
      function(match) if (length(match) >= 3) match[match_index] else NA_character_,
      character(1)
    )))
  }
  data.table(
    id = peptide_ids,
    accession = vapply(
      accession_match,
      function(match) if (length(match) >= 2) match[2] else NA_character_,
      character(1)
    ),
    start = matched_number(2),
    stop = matched_number(3)
  )
}

#' One row of the restrictions table
#'
#' @param complex,pipeline,restriction Labels of the row.
#' @param best_clustering VSClust result (`ClustOut$Bestcl`).
#' @param cannot_link The cannot-link matrix the clustering was restricted
#'   with; `NULL` for the unconstrained run.
#' @param min_form_size Integer. Members a cluster needs to count as a form.
#'
#' @return One-row data.table with the columns listed at the top of this
#'   script.
restriction_row <- function(
  complex,
  pipeline,
  restriction,
  best_clustering,
  cannot_link = NULL,
  min_form_size = 2
) {
  membership <- best_clustering$membership
  is_member <- rowMaxs(membership) >= membership_cutoff
  top_cluster <- max.col(membership, ties.method = "first")
  position <- peptide_positions(rownames(membership))
  member_top <- matrix(FALSE, nrow(membership), ncol(membership))
  member_top[cbind(seq_len(nrow(membership)), top_cluster)] <- TRUE
  count_forms <- function(member) {
    potential_forms(member, position$start, position$stop, accession = position$accession)
  }
  data.table(
    complex = complex, pipeline = pipeline, restriction = restriction,
    k = ncol(membership),
    cl_pairs = if (is.null(cannot_link)) NA_real_ else sum(cannot_link) / 2,
    violated = if (is.null(cannot_link)) {
      NA_real_
    } else {
      violations(best_clustering$cluster, cannot_link)
    },
    forms = sum(table(top_cluster[is_member]) >= min_form_size),
    members = mean(is_member),
    potential_position = count_forms(membership >= membership_cutoff),
    potential_position_1k = count_forms(membership >= 1 / ncol(membership)),
    potential_position_top = count_forms(member_top)
  )
}

#' Share of the quantification weight that comes from non-members
#'
#' @param membership Membership matrix, peptides x groups.
#' @param membership_threshold The threshold the pipeline quantifies with:
#'   memberships below it get no weight.
#' @return Median over groups of the share of a group's weight that comes
#'   from peptides with a membership below `membership_cutoff`.
nonmember_weight <- function(membership, membership_threshold) {
  weight <- membership
  weight[weight < membership_threshold] <- 0
  nonmember_share <- colSums(weight * (membership < membership_cutoff)) / colSums(weight)
  median(nonmember_share, na.rm = TRUE)
}

# ---- old pipeline ------------------------------------------------------------
old_pipeline <- load_old_implementation()
new_pipeline <- load_new_implementation()
old_data_list <- readRDS("ComplexoFinder/00-Data/H9_data.rds")
old_lookup <- readRDS("ComplexoFinder/00-LookupTables/ebi_cp_lookup.rds")
intensity_columns <- grep("^H9_B[0-9]+_D[0-9]+$", names(old_data_list[[1]]), value = TRUE)
condition_regex <- "^H9_B[0-9]+_(D[0-9]+)$"
n_replicates <- 3
n_conditions <- length(intensity_columns) / n_replicates

found_complexes <- quietly(
  old_pipeline$find_complexes(old_data_list, old_lookup, "Accession", "uniprot_id")
)
detected_proteins <- unique(sub("-1$", "", trimws(unlist(strsplit(
  unlist(lapply(old_data_list, `[[`, "Accession")), ";"
)))))
detected_proteins <- detected_proteins[!is.na(detected_proteins) & nzchar(detected_proteins)]
kept_complexes <- quietly(old_pipeline$prune_complexes_stat(
  found_complexes$complexes, old_lookup, "Accession", "uniprot_id",
  universe = detected_proteins, min_found = 2, min_total = 2, alpha = 0.05,
  adjust_method = "BH"
))$keep

old_rows <- list()
restriction_rows <- list()
for (complex in intersect(complexes, names(kept_complexes))) {
  cat("old pipeline:", complex, "\n")
  complex_data <- kept_complexes[[complex]][datatype != "NM"]
  complex_data <- old_pipeline$resolve_ambiguous_peptides(
    complex_data, accession_col = "Accession", position_col = "Position in master protein",
    lookup_table = old_lookup, id_col_lookup = "uniprot_id"
  )
  discovery_results <- quietly(old_pipeline$discover_complexoforms(
    complex_data = complex_data, intensity_cols = intensity_columns,
    group_col = "Gene name", peptide_col = "Peptide", alpha = 0.05,
    feature_adjust_method = "BH", min_total_non_na_frac = 0.40, deep_split = 2,
    min_conditions = 2, min_reps_per_condition = 2, minClusterSize = 2,
    cond_regex = condition_regex, canonical_label = "dCF0", singleton_label = "dCF-1"
  ))$peptide_results
  n_clusters <- old_pipeline$find_nclust(discovery_results)

  cannot_link <- quietly(old_pipeline$create_cannot_link_Pearson(
    data = discovery_results, complex_name = complex, id_cols = c("Gene name", "Peptide"),
    intensity_cols = intensity_columns, n_rep = n_replicates, n_cond = n_conditions,
    min_shared_cond = 10, cannot_link_th = 0.5, complex_info = old_lookup
  ))
  restrictions <- quietly(old_pipeline$vsclust_to_restrictions(
    data = discovery_results, cannotlink_matrix = cannot_link, n_clusters = n_clusters,
    n_rep = n_replicates, n_cond = n_conditions,
    id_cols = c("Gene name", "Peptide"), value_cols = intensity_columns
  ))
  constrained <- quietly(old_pipeline$vsclust_on_complex(
    data = discovery_results, id_cols = c("Gene name", "Peptide"),
    intensity_cols = intensity_columns, n_rep = n_replicates, n_cond = n_conditions,
    n_clusters = n_clusters, grouped_replicates = TRUE,
    restriction_matrix = restrictions$restriction_matrix
  ))
  initial_clustering <- restrictions$initial_clustering$ClustOut$Bestcl
  final_clustering <- constrained$ClustOut$Bestcl
  membership <- final_clustering$membership

  old_rows[[complex]] <- data.table(
    complex = complex, pipeline = "old",
    peptides = nrow(discovery_results),
    discordant = mean(discovery_results$discordant),
    k = n_clusters,
    cl_pairs = sum(cannot_link) / 2,
    violated_initial = violations(initial_clustering$cluster, cannot_link),
    violated_final = violations(final_clustering$cluster, cannot_link),
    members = mean(rowMaxs(membership) >= membership_cutoff),
    # The old pipeline quantifies with membership_threshold = 0
    nonmember_weight = if (ncol(membership) > 1) nonmember_weight(membership, 0) else 0
  )
  restriction_rows[[length(restriction_rows) + 1]] <- rbind(
    restriction_row(complex, "old", "none", initial_clustering),
    restriction_row(complex, "old", "pearson", final_clustering, cannot_link)
  )
}

# ---- new pipeline (saved results) ----------------------------------------------
#' Read the saved result of the new pipeline for one complex
read_new_result <- function(complex) {
  complex_dir <- file.path(new_results_dir, gsub("/", "_", gsub(" ", "_", complex)))
  readRDS(file.path(complex_dir, "complexoform_results.rds"))
}

new_rows <- lapply(seq_along(complexes), function(complex_index) {
  membership <- read_new_result(complexes[complex_index])$vsclust$ClustOut$Bestcl$membership
  summary_row <- new_run_summary[complex_index]
  data.table(
    complex = complexes[complex_index], pipeline = "new",
    peptides = summary_row$peptides,
    discordant = summary_row$discordant / summary_row$peptides,
    k = summary_row$n_clusters,
    cl_pairs = summary_row$cannot_link_pairs,
    violated_initial = summary_row$violated_initial,
    violated_final = summary_row$violated_final,
    members = mean(rowMaxs(membership) >= membership_cutoff),
    # The new pipeline quantifies with the members only
    nonmember_weight = if (ncol(membership) > 1) {
      nonmember_weight(membership, membership_cutoff)
    } else {
      0
    }
  )
})

realdata_metrics <- rbind(rbindlist(old_rows), rbindlist(new_rows))
fwrite(realdata_metrics, file.path(results_dir, "realdata_metrics.csv"))
cat("\nWrote", file.path(results_dir, "realdata_metrics.csv"), "\n")

# ---- new pipeline: every cannot-link rule ----------------------------------------
# With the clustering settings of NewImplementation/RunComplexoFinder.R. The
# pipeline's own rule (groups, default tolerance) is rebuilt here too and
# checked against the saved result.
discovery <- readRDS(file.path(new_results_dir, "discovery.rds"))
peptide_key <- function(peptide_table) {
  paste(
    peptide_table$complex_name, peptide_table[["Gene name"]], peptide_table$Peptide,
    sep = "|"
  )
}
discovery_keys <- peptide_key(discovery$peptide_results)
id_columns <- c("Gene name", "Peptide")
# "tolerance" and "groups" with the pipeline's default tolerance
rules <- list(
  pearson = list(method = "pearson"),
  calibrated = list(method = "calibrated"),
  zrmsd = list(method = "zrmsd"),
  ccc = list(method = "ccc"),
  tolerance = list(method = "tolerance"),
  groups = list(method = "groups"),
  "groups, no tolerance" = list(method = "groups", tolerance = 0)
)

for (complex in complexes) {
  cat("new pipeline, cannot-link rules:", complex, "\n")
  saved_result <- read_new_result(complex)
  if (is.null(saved_result$constraints)) {
    # One cluster, or the run was made without constraints
    restriction_rows[[length(restriction_rows) + 1]] <-
      restriction_row(complex, "new", "none", saved_result$vsclust$ClustOut$Bestcl)
    next
  }
  complex_data <- saved_result$dCF
  initial_clustering <- saved_result$constraints$initial_clustering$ClustOut$Bestcl
  initial_centers <- initial_clustering$centers
  rownames(initial_centers) <- seq_len(nrow(initial_centers))
  discovery_rows <- match(peptide_key(complex_data), discovery_keys)
  stopifnot(!anyNA(discovery_rows))
  noise_model <- list(
    squared_standard_error = discovery$model$se2[discovery_rows, , drop = FALSE],
    between_peptide_variance = discovery$tau2
  )

  restriction_rows[[length(restriction_rows) + 1]] <-
    restriction_row(complex, "new", "none", initial_clustering)
  for (rule in names(rules)) {
    cannot_link <- quietly(do.call(new_pipeline$create_cannot_link, c(
      list(
        data = complex_data, id_columns = id_columns,
        intensity_columns = intensity_columns, condition_regex = condition_regex,
        min_shared_conditions = 10, noise_model = noise_model
      ),
      rules[[rule]]
    )))
    best_clustering <- if (sum(cannot_link) == 0) {
      initial_clustering
    } else {
      quietly(new_pipeline$vsclust_on_complex(
        data = complex_data, id_columns = id_columns,
        intensity_columns = intensity_columns, condition_regex = condition_regex,
        n_clusters = nrow(initial_centers), scaling = "center", seed = 1,
        restriction_matrix = restrictions_from_clustering(initial_clustering, cannot_link),
        initial_centers = initial_centers
      ))$ClustOut$Bestcl
    }
    # "groups" with the default tolerance is the pipeline's own rule
    if (rule == "groups") {
      saved_membership <- saved_result$vsclust$ClustOut$Bestcl$membership
      rebuilt_membership <- best_clustering$membership[rownames(saved_membership), ]
      if (!isTRUE(all.equal(rebuilt_membership, saved_membership, check.attributes = FALSE))) {
        warning(
          complex, ": the clustering rebuilt here with the pipeline's rule differs from the ",
          "saved one. The settings of the saved run may differ from those assumed in this script."
        )
      }
    }
    restriction_rows[[length(restriction_rows) + 1]] <-
      restriction_row(complex, "new", rule, best_clustering, cannot_link)
  }
}

fwrite(rbindlist(restriction_rows), file.path(results_dir, "realdata_restrictions.csv"))
cat("Wrote", file.path(results_dir, "realdata_restrictions.csv"), "\n")
