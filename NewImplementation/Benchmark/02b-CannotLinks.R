# ============================================================================
# 02b-CannotLinks.R
# ============================================================================
# The cannot-link rules themselves, on ProteoMaker, without any clustering.
#
# A cannot-link says: these two overlapping peptides are different forms of the
# same stretch of the protein. In ProteoMaker that can be checked pair by pair:
# two peptides truly differ when their expected profiles differ, i.e. when
# they are in different truth groups.
#
# Every rule scores the same candidate pairs: two peptides of one protein
# that overlap by more than half of the shorter one and share at least
# `min_shared_conditions` conditions. The peptides are those the pipeline
# would cluster (tested, seen in >= `min_conditions_clustering` conditions),
# in every protein of every dataset. The noise model comes from the new
# discovery, as in the pipeline.
#
#   pearson      Pearson r of the condition means < 0.5
#   calibrated   pair by pair: profiles differ beyond the noise model
#   zrmsd        root mean square difference of robust z-scores > 1
#   ccc          Lin's concordance < 0.5
#   tolerance 0.25, tolerance 0.5
#                pair by pair: profiles differ by more than that many log2
#                beyond the noise model. (The pipeline's default tolerance,
#                the square root of the between-peptide variance, is 0 on
#                ProteoMaker, so fixed tolerances are used here.)
#   groups, groups 0.25, groups 0.5
#                pattern groups per position, without and with a tolerance
#
# Per dataset and rule:
#   candidates    candidate pairs
#   different     share of them that truly differ
#   linked        share linked
#   precision     of the linked pairs, the share that truly differ
#   recall        of the truly different pairs, the share linked
#   false_links   of the pairs that do not differ, the share linked
#   AUROC         how well the rule's score ranks the truly different pairs
#                 above the others, whatever the threshold
#   precision_at_n  precision of each rule's n most confident pairs, n being
#                 the number of links of the rule that makes fewest in the
#                 dataset: the rules at the same number of links
#
# The dataset "all" pools the six datasets.
#
# Writes to results/:
#   cannotlink_pairs.csv                  one row per candidate pair and rule
#   cannotlink_metrics.csv                the metrics above
#   cannotlink_linked_by_difference.csv   share of pairs linked, by the size of
#                                         the true difference between the two
#                                         expected profiles (RMS, log2)
#
# Run from the project root (~3 min):
#   Rscript NewImplementation/Benchmark/02b-CannotLinks.R
# ============================================================================

suppressMessages(library(data.table))
setwd(here::here())
source("NewImplementation/Benchmark/BenchmarkData.R")
source("NewImplementation/Benchmark/BenchmarkFunctions.R")

results_dir <- "NewImplementation/Benchmark/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

new_pipeline <- load_new_implementation()

# As in 02-Clustering.R
min_conditions_clustering <- 4
min_shared_conditions <- 3
# Rule label, create_cannot_link() method, tolerance (log2; NA: the default)
rules <- data.table(
  rule = c("pearson", "calibrated", "zrmsd", "ccc", "tolerance 0.25", "tolerance 0.5",
           "groups", "groups 0.25", "groups 0.5"),
  method = c("pearson", "calibrated", "zrmsd", "ccc", "tolerance", "tolerance",
             "groups", "groups", "groups"),
  tolerance_log2 = c(NA, NA, NA, NA, 0.25, 0.5, 0, 0.25, 0.5)
)
# TRUE: a higher score means "more different"
higher_score_is_different <- setNames(rules$method == "zrmsd", rules$rule)

pairs <- list()
for (dataset_name in proteomaker_datasets) {
  dataset <- load_proteomaker_dataset(dataset_name)
  discovery <- run_discovery_new(new_pipeline, dataset)$discovery
  id_columns <- c(dataset$assembly_column, dataset$peptide_column)
  all_peptide_ids <- paste(
    dataset$data[[dataset$assembly_column]],
    dataset$data[[dataset$peptide_column]],
    sep = "_"
  )
  # Expected profiles, each centered on its own mean: what the rules compare
  centered_expected <- dataset$expected - rowMeans(dataset$expected)
  conditions_observed <- rowSums(is.finite(discovery$model$cond_mean))
  is_clustered <- discovery$model$tested &
    conditions_observed >= min_conditions_clustering

  rows_by_protein <- split(which(is_clustered), dataset$data$Accession[is_clustered])
  rows_by_protein <- rows_by_protein[lengths(rows_by_protein) >= 2]
  cat(sprintf(
    "[%s] %d proteins with at least two clustered peptides, between-peptide sd = %.3f\n",
    dataset_name, length(rows_by_protein), sqrt(discovery$tau2)
  ))

  for (protein in names(rows_by_protein)) {
    protein_rows <- rows_by_protein[[protein]]
    protein_data <- discovery$peptide_results[protein_rows]
    noise_model <- list(
      squared_standard_error = discovery$model$se2[protein_rows, , drop = FALSE],
      between_peptide_variance = discovery$tau2
    )
    for (rule_index in seq_len(nrow(rules))) {
      tolerance <- if (is.na(rules$tolerance_log2[rule_index])) {
        NULL
      } else {
        rules$tolerance_log2[rule_index]
      }
      scored_pairs <- attr(quietly(new_pipeline$create_cannot_link(
        protein_data,
        id_columns = id_columns,
        intensity_columns = dataset$intensity_columns,
        condition_regex = dataset$condition_regex,
        method = rules$method[rule_index],
        tolerance = tolerance,
        min_shared_conditions = min_shared_conditions,
        accession_column = dataset$assembly_column,
        position_column = "Position",
        noise_model = noise_model
      )), "pairs")
      if (is.null(scored_pairs) || nrow(scored_pairs) == 0) next
      first_row <- match(scored_pairs$peptide_i, all_peptide_ids)
      second_row <- match(scored_pairs$peptide_j, all_peptide_ids)
      pairs[[length(pairs) + 1]] <- data.table(
        dataset = dataset_name, protein = protein, method = rules$rule[rule_index],
        peptide_i = scored_pairs$peptide_i, peptide_j = scored_pairs$peptide_j,
        shared = scored_pairs$shared, score = scored_pairs$score,
        linked = scored_pairs$linked,
        different = dataset$truth$truth_group[first_row] !=
          dataset$truth$truth_group[second_row],
        # Size of the true difference: root mean square over conditions (log2)
        true_distance = sqrt(rowMeans(
          (centered_expected[first_row, , drop = FALSE] -
            centered_expected[second_row, , drop = FALSE])^2
        ))
      )
    }
  }
}
pairs <- rbindlist(pairs)
fwrite(pairs, file.path(results_dir, "cannotlink_pairs.csv"))

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
# A higher confidence means "more different", whatever the rule's score
pairs[, confidence := ifelse(higher_score_is_different[method], score, -score)]
# The rules at the same number of links: that of the rule making fewest
pairs[, n_matched := min(tapply(linked, method, sum)), by = dataset]

#' Metrics of one rule on one set of candidate pairs
#'
#' @param rule_pairs data.table of candidate pairs with `linked`, `different`,
#'   `confidence` and `n_matched`.
#' @return One-row data.table with the metrics listed at the top of this
#'   script.
pair_metrics <- function(rule_pairs) {
  ranked <- rule_pairs[is.finite(confidence)][order(-confidence)]
  n_matched <- rule_pairs$n_matched[1]
  data.table(
    candidates = nrow(rule_pairs),
    different = mean(rule_pairs$different),
    linked = mean(rule_pairs$linked),
    precision = if (any(rule_pairs$linked)) {
      mean(rule_pairs$different[rule_pairs$linked])
    } else {
      NA_real_
    },
    recall = if (any(rule_pairs$different)) {
      mean(rule_pairs$linked[rule_pairs$different])
    } else {
      NA_real_
    },
    false_links = if (any(!rule_pairs$different)) {
      mean(rule_pairs$linked[!rule_pairs$different])
    } else {
      NA_real_
    },
    AUROC = auc_rank(-ranked$confidence, ranked$different),
    n_matched = n_matched,
    precision_at_n = if (n_matched > 0 && nrow(ranked) >= n_matched) {
      mean(ranked$different[seq_len(n_matched)])
    } else {
      NA_real_
    }
  )
}
pooled_pairs <- copy(pairs)[, dataset := "all"]
pooled_pairs[, n_matched := min(tapply(linked, method, sum))]
metrics <- rbind(
  pairs[, pair_metrics(.SD), by = .(dataset, method)],
  pooled_pairs[, pair_metrics(.SD), by = .(dataset, method)]
)
metrics[, method := factor(method, levels = rules$rule)]
setorder(metrics, dataset, method)
fwrite(metrics, file.path(results_dir, "cannotlink_metrics.csv"))

# Share linked by the size of the true difference
difference_levels <- c("none", "below 0.25", "0.25 to 0.5", "0.5 to 1", "above 1")
pairs[, difference := as.character(cut(
  true_distance, c(-Inf, 0.25, 0.5, 1, Inf),
  labels = difference_levels[-1]
))]
pairs[different == FALSE, difference := "none"]
pairs[, difference := factor(difference, levels = difference_levels)]
pairs[, method := factor(method, levels = rules$rule)]
linked_by_difference <- dcast(
  pairs[, .(pairs = .N, linked = mean(linked)), by = .(difference, method)],
  difference + pairs ~ method, value.var = "linked"
)
fwrite(linked_by_difference, file.path(results_dir, "cannotlink_linked_by_difference.csv"))

print(metrics[, lapply(.SD, function(column) if (is.double(column)) round(column, 3) else column)])
cat("\nWrote", file.path(results_dir, "cannotlink_metrics.csv"), "\n")
