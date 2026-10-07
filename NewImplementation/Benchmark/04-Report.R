# ============================================================================
# 04-Report.R
# ============================================================================
# Turns the outputs of 01-03 (and 02b) into the comparison of the old and the
# new implementation:
#
#   results/REPORT.md                        the tables below, in one document
#   results/summary_discovery.csv            old vs new per dataset
#   results/summary_clustering.csv           every configuration per dataset
#   results/summary_clustering_paired.csv    new (pipeline) vs old (pipeline),
#                                            protein by protein
#   results/summary_forms.csv                number of forms per configuration
#   results/summary_realdata.csv             H9 complexes, per pipeline
#   results/figure_discovery.png
#   results/figure_clustering.png
#
# Uses whichever of the three result files exist.
#
# Run from the project root:
#   Rscript NewImplementation/Benchmark/04-Report.R
# ============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
})
setwd(here::here())

results_dir <- "NewImplementation/Benchmark/results"
#' Read a results table if it exists
#'
#' @param file_name Name of a csv file in the results directory.
#' @return data.table, or `NULL` if there is no such file.
read_if <- function(file_name) {
  file_path <- file.path(results_dir, file_name)
  if (file.exists(file_path)) fread(file_path) else NULL
}
discovery <- read_if("discovery_metrics.csv")
clustering <- read_if("clustering_metrics.csv")
realdata <- read_if("realdata_metrics.csv")
restrictions <- read_if("realdata_restrictions.csv")
link_metrics <- read_if("cannotlink_metrics.csv")
link_by_difference <- read_if("cannotlink_linked_by_difference.csv")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Format a table as a markdown table
#'
#' @param table_data data.table or data.frame.
#' @param digits Decimals shown for non-integer numbers.
#' @return Character vector of markdown lines; `NA` is shown as a dash.
md_table <- function(table_data, digits = 2) {
  table_data <- as.data.table(table_data)
  format_column <- function(column) {
    if (is.numeric(column) && !is.integer(column)) {
      ifelse(is.na(column), "–", formatC(column, format = "f", digits = digits))
    } else {
      ifelse(is.na(column), "–", as.character(column))
    }
  }
  table_rows <- do.call(paste, c(lapply(table_data, format_column), sep = " | "))
  c(
    paste0("| ", paste(names(table_data), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(table_data)), collapse = "|"), "|"),
    paste0("| ", table_rows, " |")
  )
}

#' Mean of the non-missing values
#'
#' @param values Numeric vector.
#' @return The mean without `NA`s, or `NA` when every value is missing.
mean_na <- function(values) {
  if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
}

#' Sign test of one set of values against another, as a table cell
#'
#' @param first_values,second_values Numeric vectors, paired.
#' @param lower_is_better Logical. Whether a lower value is the better one.
#' @return Character: in how many pairs the first is better / worse, and the
#'   sign-test p-value in brackets.
sign_cell <- function(first_values, second_values, lower_is_better = FALSE) {
  both_finite <- is.finite(first_values) & is.finite(second_values)
  advantage <- if (lower_is_better) {
    second_values[both_finite] - first_values[both_finite]
  } else {
    first_values[both_finite] - second_values[both_finite]
  }
  better <- sum(advantage > 1e-9)
  worse <- sum(advantage < -1e-9)
  p_value <- if (better + worse > 0) {
    stats::binom.test(better, better + worse)$p.value
  } else {
    NA_real_
  }
  p_text <- if (is.na(p_value)) {
    "–"
  } else if (p_value < 0.001) {
    "<0.001"
  } else {
    formatC(p_value, format = "f", digits = 3)
  }
  sprintf("%d / %d (%s)", better, worse, p_text)
}

#' Log10 of the positive values
#'
#' @param values Vector that can be read as numbers.
#' @return log10 of the values; `NA` for values that are not positive.
log_positive <- function(values) {
  values <- as.numeric(values)
  values[!is.na(values) & values <= 0] <- NA_real_
  log10(values)
}

# Colours: new = blue, old = orange (a colour-blind-safe pair), published
# tools = grey.
new_colour <- "#2a78d6"
old_colour <- "#eb6834"
reference_colour <- "#898781"

#' ggplot theme of the benchmark figures
theme_bench <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "#fcfcfb", colour = NA),
      panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(colour = "#52514e"),
      axis.title = element_blank(),
      strip.text = element_text(colour = "#0b0b0b", face = "bold", hjust = 0),
      strip.text.y = element_text(angle = 0),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      plot.title = element_text(colour = "#0b0b0b", face = "bold"),
      plot.subtitle = element_text(colour = "#52514e"),
      panel.spacing = unit(1.2, "lines")
    )
}

report <- c(
  "# Benchmark report: old vs new implementation",
  "",
  paste0("*Generated ", Sys.Date(), " by `NewImplementation/Benchmark/04-Report.R`.*"),
  "",
  "- **old** = the pipeline in `ComplexoFinder/`, run as its `RunComplexoFinder.R` runs it.",
  "- **new** = the pipeline in `NewImplementation/`.",
  "",
  "Both are run by the same harness on the same data and scored with the same metrics.",
  "What each benchmark and metric means is described in [../README.md](../README.md).",
  ""
)

# ---------------------------------------------------------------------------
# 1. Discovery
# ---------------------------------------------------------------------------
if (!is.null(discovery)) {
  dataset_levels <- unique(discovery$dataset)
  discovery_summary <- discovery[, .(
    benchmark, dataset, method, coverage, AUROC = auroc,
    within, false_alarms = fpr, power = tpr, false_calls = fdp
  )]
  fwrite(discovery_summary, file.path(results_dir, "summary_discovery.csv"))

  report <- c(
    report,
    "## 1. Discovery of discordant peptides",
    "",
    "- **coverage**: share of peptides the method scores at all.",
    "- **within**: AUROC inside proteins, i.e. can it tell a discordant peptide from its own siblings.",
    "- **false_alarms**: share of peptides called in proteins with no discordant peptide.",
    "- **power**: share of discordant peptides called.",
    "- **false_calls**: share of the called peptides that are not discordant. This is what the nominal 5% is meant to bound.",
    "",
    "`new_global_bh` is the new implementation with Benjamini-Hochberg across all peptides only, without the Bonferroni step within each protein.",
    "",
    "PeCorA, ProteoForge and COPF are saved results of the published tools, scored the same way, for reference.",
    ""
  )
  for (benchmark_name in unique(discovery_summary$benchmark)) {
    report <- c(
      report,
      paste0("### ", benchmark_name), "",
      md_table(discovery_summary[benchmark == benchmark_name, !"benchmark"], 3), ""
    )
  }

  # Figure: old -> new per dataset, one panel per metric
  metrics_long <- melt(
    discovery_summary[, .(
      benchmark, dataset, method,
      `Within-protein AUROC\n(higher is better)` = within,
      `False alarms\n(lower is better)` = false_alarms,
      `Power\n(higher is better)` = power
    )],
    id.vars = c("benchmark", "dataset", "method"), variable.name = "metric"
  )
  metrics_long[, dataset := factor(dataset, levels = rev(dataset_levels))]
  metrics_long[, benchmark := factor(benchmark, levels = unique(discovery_summary$benchmark))]
  old_and_new <- dcast(
    metrics_long[method %in% c("old", "new")],
    benchmark + dataset + metric ~ method
  )
  figure <- ggplot() +
    geom_segment(data = old_and_new, aes(x = old, xend = new, y = dataset, yend = dataset),
                 colour = "#c3c2b7", linewidth = 0.7) +
    geom_point(data = metrics_long[!method %in% c("old", "new", "new_global_bh")],
               aes(x = value, y = dataset, shape = "Published tools (saved results)"),
               colour = reference_colour, size = 2.2, stroke = 0.7, na.rm = TRUE) +
    geom_point(data = metrics_long[method %in% c("old", "new")],
               aes(x = value, y = dataset, colour = method), size = 3, na.rm = TRUE) +
    geom_point(data = metrics_long[method == "new_global_bh"],
               aes(x = value, y = dataset, shape = "New, global BH only"),
               colour = new_colour, fill = "#fcfcfb", size = 2.6, stroke = 1, na.rm = TRUE) +
    scale_colour_manual(values = c(old = old_colour, new = new_colour),
                        breaks = c("old", "new"),
                        labels = c(old = "Old implementation", new = "New implementation")) +
    scale_shape_manual(values = c("Published tools (saved results)" = 1,
                                  "New, global BH only" = 21)) +
    guides(shape = guide_legend(override.aes = list(
      colour = c(new_colour, reference_colour), fill = "#fcfcfb"
    ))) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    facet_grid(benchmark ~ metric, scales = "free_y", space = "free_y") +
    labs(title = "Discovery of discordant peptides",
         subtitle = "Each line joins the old and the new implementation on one dataset") +
    theme_bench()
  ggsave(file.path(results_dir, "figure_discovery.png"), figure, width = 9, height = 5.2, dpi = 200)
  report <- c(report, "![Discovery](figure_discovery.png)", "")
}

# ---------------------------------------------------------------------------
# 2. Clustering and quantification (ProteoMaker)
# ---------------------------------------------------------------------------
if (!is.null(clustering)) {
  clustering[, k_error := abs(k - k_true)]
  # ARI is only informative when there is more than one truth group
  clustering[, ari_structured := ifelse(truth_groups > 1, ari, NA_real_)]
  # Forms: clusters with members against truth groups, both of >= 2 peptides
  clustering[, forms_error := abs(forms_found - k_true)]
  # Potential proteoforms (old clustering benchmark): compared on the log10
  # scale as there, where a clustering without members (0) drops out
  for (column in grep("^potential_", names(clustering), value = TRUE)) {
    set(clustering, j = paste0("log_", column), value = log_positive(clustering[[column]]))
  }
  dataset_levels <- unique(clustering$dataset)

  #' Share of cannot-link pairs in one cluster, pooled over proteins
  pooled_violations <- function(violated_share, n_pairs) {
    if (sum(n_pairs, na.rm = TRUE) > 0) {
      weighted.mean(violated_share, n_pairs, na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  #' Summary of the clustering metrics per group of rows
  #'
  #' @param metrics data.table of clustering metrics, one row per protein and
  #'   configuration.
  #' @param by Character vector of columns to summarise by.
  summarise <- function(metrics, by) {
    metrics[, .(
      proteins = .N,
      k = mean(k), k_true = mean(k_true), k_error = mean(k_error),
      assigned = mean(assigned),
      pair_precision = mean_na(pair_precision),
      pair_recall = mean_na(pair_recall),
      ARI = mean_na(ari_structured),
      profile_R2 = mean_na(profile_r2),
      top_fraction = mean_na(top_fraction),
      potential = mean_na(log_potential_proteoforms),
      potential_1k = mean_na(log_potential_proteoforms_1k),
      potential_top = mean_na(log_potential_proteoforms_top),
      pattern = mean_na(log_potential_pattern),
      pattern_1k = mean_na(log_potential_pattern_1k),
      pattern_top = mean_na(log_potential_pattern_top),
      pattern_tol_1k = mean_na(log_potential_pattern_tol_1k),
      pattern_tol_top = mean_na(log_potential_pattern_tol_top),
      position = mean_na(log_potential_position),
      position_1k = mean_na(log_potential_position_1k),
      position_top = mean_na(log_potential_position_top),
      sum_unique = mean_na(as.numeric(sum_unique_proteoforms)),
      sum_unique_1k = mean_na(as.numeric(sum_unique_proteoforms_1k)),
      forms_present = sum(k_true), forms_found = sum(forms_found),
      forms_recovered = sum(forms_recovered), forms_error = mean(forms_error),
      cl_pairs = sum(cl_pairs, na.rm = TRUE),
      cl_violated = pooled_violations(cl_violated, cl_pairs),
      quant_all_plain = mean_na(quant_all_plain),
      quant_all_centred = mean_na(quant_all_centred),
      quant_members_plain = mean_na(quant_members_plain),
      quant_members_centred = mean_na(quant_members_centred)
    ), by = by]
  }
  summary_by_dataset <- summarise(clustering, c("dataset", "config"))
  summary_overall <- summarise(clustering, "config")
  fwrite(summary_by_dataset, file.path(results_dir, "summary_clustering.csv"))

  # Protein by protein: new (pipeline) against old (pipeline)
  paired_metrics <- c("k_error", "assigned", "pair_precision", "pair_recall",
                      "ari_structured", "profile_r2", "top_fraction",
                      "log_potential_proteoforms", "log_potential_proteoforms_1k",
                      "log_potential_proteoforms_top",
                      "log_potential_position", "log_potential_position_1k",
                      "log_potential_position_top",
                      "forms_error", "forms_recovered")
  lower_is_better_metrics <- c(
    "k_error", "forms_error", grep("^log_potential", paired_metrics, value = TRUE)
  )
  pipelines_wide <- dcast(
    clustering[config %in% c("old (pipeline)", "new (pipeline)")],
    dataset + protein ~ config, value.var = paired_metrics
  )
  paired_summary <- rbindlist(lapply(paired_metrics, function(metric) {
    old_values <- pipelines_wide[[paste0(metric, "_old (pipeline)")]]
    new_values <- pipelines_wide[[paste0(metric, "_new (pipeline)")]]
    both_finite <- is.finite(old_values) & is.finite(new_values)
    # Positive where the new pipeline is better
    improvement <- if (metric %in% lower_is_better_metrics) {
      old_values[both_finite] - new_values[both_finite]
    } else {
      new_values[both_finite] - old_values[both_finite]
    }
    differs <- abs(improvement) > 1e-9
    data.table(
      metric = metric, proteins = sum(both_finite),
      old = mean(old_values[both_finite]), new = mean(new_values[both_finite]),
      new_better = sum(improvement > 1e-9), same = sum(!differs),
      new_worse = sum(improvement < -1e-9),
      p_sign = if (any(differs)) {
        stats::binom.test(sum(improvement > 1e-9), sum(differs))$p.value
      } else {
        NA_real_
      }
    )
  }))
  fwrite(paired_summary, file.path(results_dir, "summary_clustering_paired.csv"))

  cluster_columns <- c("assigned", "pair_precision", "pair_recall", "ARI", "profile_R2",
                       "top_fraction", "cl_violated")

  # Number of forms under each way of restricting the clustering, over the
  # proteins that have more than one truth group (as for the ARI)
  restriction_configs <- c("old, unconstrained", "old (pipeline)", "new, unconstrained",
                           "new (pipeline)", "new, calibrated links", "new, ZRMSD links",
                           "new, CCC links", "new, tolerance links (0.25)", "new, group links",
                           "new, group links (0.25)")
  restriction_configs <- intersect(restriction_configs, unique(clustering$config))
  structured <- clustering[truth_groups > 1 & config %in% restriction_configs]
  structured[, config := factor(config, levels = restriction_configs)]
  forms_summary <- structured[, .(
    proteins = .N,
    cl_pairs = sum(cl_pairs, na.rm = TRUE),
    present = sum(k_true), requested = sum(k), found = sum(forms_found),
    recovered = sum(forms_recovered),
    recovered_share = sum(forms_recovered) / sum(k_true),
    found_right = sum(forms_recovered) / sum(forms_found),
    count_error = mean(forms_error),
    profile_R2 = mean_na(profile_r2)
  ), keyby = config]
  fwrite(forms_summary, file.path(results_dir, "summary_forms.csv"))
  forms_by_dataset <- dcast(
    structured[,
      .(present_found_recovered = sprintf(
        "%d / %d / %d", sum(k_true), sum(forms_found), sum(forms_recovered)
      )),
      by = .(dataset, config)
    ],
    dataset ~ config, value.var = "present_found_recovered"
  )[match(intersect(dataset_levels, structured$dataset), dataset)]

  # One restriction against another, protein by protein (all proteins)
  by_protein <- dcast(clustering[config %in% restriction_configs], dataset + protein ~ config,
                      value.var = c("forms_recovered", "forms_error", "profile_r2",
                                    "log_potential_proteoforms", "log_potential_proteoforms_1k",
                                    "log_potential_proteoforms_top",
                                    "log_potential_position", "log_potential_position_1k",
                                    "log_potential_position_top",
                                    "log_potential_pattern_1k", "log_potential_pattern_top",
                                    "log_potential_pattern_tol_1k", "log_potential_pattern_tol_top"))
  comparisons <- list(
    c("new (pipeline)", "new, unconstrained"), c("new, calibrated links", "new, unconstrained"),
    c("new, ZRMSD links", "new, unconstrained"), c("new, CCC links", "new, unconstrained"),
    c("new, tolerance links (0.25)", "new, unconstrained"), c("new, group links", "new, unconstrained"),
    c("new, group links (0.25)", "new, unconstrained"),
    c("new, calibrated links", "new (pipeline)"), c("new, ZRMSD links", "new (pipeline)"),
    c("new, CCC links", "new (pipeline)"), c("new, tolerance links (0.25)", "new (pipeline)"),
    c("new, group links", "new (pipeline)"), c("new, group links (0.25)", "new (pipeline)"),
    c("new, tolerance links (0.25)", "new, calibrated links"),
    c("new, group links", "new, calibrated links"),
    c("new, group links (0.25)", "new, tolerance links (0.25)"),
    c("old (pipeline)", "old, unconstrained"), c("new (pipeline)", "old (pipeline)")
  )
  comparisons <- Filter(
    function(comparison) all(comparison %in% restriction_configs),
    comparisons
  )
  #' Sign-test cell for one metric: the first configuration of a comparison
  #' against the second
  compare <- function(comparison, metric, lower_is_better = FALSE) {
    sign_cell(
      by_protein[[paste0(metric, "_", comparison[1])]],
      by_protein[[paste0(metric, "_", comparison[2])]],
      lower_is_better = lower_is_better
    )
  }
  forms_paired <- rbindlist(lapply(comparisons, function(comparison) {
    data.table(
      configuration = comparison[1], against = comparison[2],
      forms_recovered = compare(comparison, "forms_recovered"),
      count_error = compare(comparison, "forms_error", lower_is_better = TRUE),
      profile_R2 = compare(comparison, "profile_r2")
    )
  }))
  potential_paired <- rbindlist(lapply(comparisons, function(comparison) {
    compare_potential <- function(metric) {
      compare(comparison, metric, lower_is_better = TRUE)
    }
    data.table(
      configuration = comparison[1], against = comparison[2],
      potential = compare_potential("log_potential_proteoforms"),
      potential_1k = compare_potential("log_potential_proteoforms_1k"),
      potential_top = compare_potential("log_potential_proteoforms_top"),
      position = compare_potential("log_potential_position"),
      position_1k = compare_potential("log_potential_position_1k"),
      position_top = compare_potential("log_potential_position_top"),
      pattern_1k = compare_potential("log_potential_pattern_1k"),
      pattern_top = compare_potential("log_potential_pattern_top"),
      pattern_tol_1k = compare_potential("log_potential_pattern_tol_1k"),
      pattern_tol_top = compare_potential("log_potential_pattern_tol_top"),
      profile_R2 = compare(comparison, "profile_r2")
    )
  }))
  potential_columns <- c("config", "assigned", "potential", "potential_1k", "potential_top",
                         "position", "position_1k", "position_top", "sum_unique", "sum_unique_1k")
  pattern_columns <- c("config", "pattern_1k", "pattern_top", "pattern_tol_1k",
                       "pattern_tol_top", "profile_R2")
  potential_comparison_columns <- c("configuration", "against", "potential", "potential_1k",
                                    "potential_top", "position", "position_1k", "position_top")
  pattern_comparison_columns <- c("configuration", "against", "pattern_1k", "pattern_top",
                                  "pattern_tol_1k", "pattern_tol_top", "profile_R2")
  summary_with_missing <- summarise(
    clustering[grepl("NA$", dataset) & !grepl("NoNA$", dataset)],
    "config"
  )
  pipeline_summary <- summary_by_dataset[
    config %in% c("old (pipeline)", "new (pipeline)"),
    c("dataset", "config", "k", "k_true", "k_error", cluster_columns),
    with = FALSE
  ]

  # The cannot-link rules pair by pair (02b-CannotLinks.R)
  link_section <- character()
  if (!is.null(link_metrics)) {
    link_columns <- c("method", "candidates", "different", "linked", "precision", "recall",
                      "false_links", "AUROC", "precision_at_n")
    links_by_dataset <- dcast(
      link_metrics[
        dataset != "all",
        .(dataset, method,
          false_links_and_auroc = sprintf("%.2f / %.2f", false_links, AUROC))
      ],
      dataset ~ method, value.var = "false_links_and_auroc"
    )[match(intersect(dataset_levels, link_metrics$dataset), dataset)]
    link_section <- c(
      "### The cannot-link rules, pair by pair",
      "",
      "Without any clustering: every pair of overlapping peptides of a protein that a rule can link, in every protein of the six datasets, against whether the two truly differ (different expected profiles).",
      "",
      "- **candidates**: pairs a rule can link. **different**: share of them that truly differ.",
      "- **linked**: share the rule links. **precision**: of the linked pairs, the share that truly differ. **recall**: of the truly different pairs, the share linked.",
      "- **false_links**: of the pairs that do not differ, the share linked.",
      "- **AUROC**: how well the rule's score ranks the truly different pairs above the others, whatever the threshold (0.5 = no better than chance).",
      "- **precision_at_n**: precision of each rule's n most confident pairs, n being the number of links of the rule that makes fewest.",
      "",
      md_table(link_metrics[dataset == "all", ..link_columns], 3), "",
      "Per dataset, as false_links / AUROC:",
      "",
      md_table(links_by_dataset), ""
    )
    if (!is.null(link_by_difference)) {
      link_section <- c(
        link_section,
        "Share of pairs linked, by the size of the true difference between the two expected profiles (RMS over conditions, log2):",
        "",
        md_table(link_by_difference), ""
      )
    }
  }

  report <- c(
    report,
    "## 2. Clustering (ProteoMaker, end to end)",
    "",
    paste0("The ", uniqueN(clustering$protein), " proteins with most peptides, in each of ",
           uniqueN(clustering$dataset), " datasets. Values are means over proteins."),
    "",
    "- **k_error**: |clusters used − truth groups with at least two peptides|.",
    "- **assigned**: share of a protein's peptides that are a member of a cluster (membership ≥ 0.5).",
    "- **pair_precision**: of the peptide pairs put in one cluster, the share that belong to one truth group. Low means clusters mix groups.",
    "- **pair_recall**: of the peptide pairs in one truth group, the share put in one cluster. Low means groups are split or left unassigned.",
    "- **ARI**: adjusted Rand index against the truth groups, over proteins with more than one truth group; an unassigned peptide is a group of its own.",
    "- **profile_R2**: share of the true profile differences that the clusters explain.",
    "- **top_fraction**: the old benchmark's headline: of the peptides carrying proteoform 1, the largest share in one cluster.",
    "- **cl_violated**: share of cannot-link pairs that sit in one cluster.",
    "",
    "### Pipeline against pipeline",
    "",
    md_table(pipeline_summary), "",
    "### Protein by protein, all datasets",
    "",
    "`p_sign` is a sign test on the proteins where the two differ.",
    "",
    md_table(paired_summary, 3), "",
    "### Does constraining help?",
    "",
    md_table(summary_overall[
      match(restriction_configs, config), c("config", cluster_columns), with = FALSE
    ]), "",
    link_section,
    "### Potential proteoforms (the measure of the old clustering benchmark)",
    "",
    "`Benchmark/03-Clustering` reads a cluster as one form and counts how many forms it leaves open. The cluster's members are laid out along the sequence and merged into stretches of overlapping peptides; the cluster's potential forms are the product over those stretches of the number of alternatives in each, and a clustering's value is the sum over its clusters. The functions here give the same values as the old benchmark's on the same clusterings.",
    "",
    "- **potential**: the old *Total potential proteoforms*. The alternatives in a stretch are the distinct true proteoform IDs its member peptides carry.",
    "- **position**: the old *Position-only potential proteoforms*. The alternatives are the member peptides of the stretch; this needs no truth.",
    "- **sum_unique**: the old *Sum of unique proteoforms*: distinct proteoform IDs among a cluster's members, summed over clusters.",
    "- **_1k**: the old benchmark's flexible cutoff. A peptide is a member of every cluster it has a membership of at least 1/k in, k being the number of clusters. Every clustered peptide is then a member of at least one cluster and can be a member of several.",
    "- **potential_top**, **position_top**: every clustered peptide counted once, in the cluster it has the highest membership in.",
    "",
    "All are means over proteins, **lower is better**, and potential and position are on the log10 scale, as in the old benchmark. With the fixed cutoff only members (membership ≥ 0.5) count, so a clustering that assigns fewer peptides scores better and a cluster without members counts nothing; read those columns next to **assigned**. The `_1k` and `_top` columns count every clustered peptide.",
    "",
    md_table(summary_overall[
      match(restriction_configs, config), potential_columns, with = FALSE
    ]), "",
    "The datasets with missing values only (LowNA, MedNA, HighNA), which is what the old benchmark was run on:",
    "",
    md_table(summary_with_missing[
      match(restriction_configs, config), potential_columns, with = FALSE
    ]), "",
    "Protein by protein, over all proteins: in how many the first configuration is better (lower; higher for profile_R2) / worse than the second, with the sign-test p-value in brackets.",
    "",
    md_table(potential_paired[, ..potential_comparison_columns]), "",
    "### Potential forms by true pattern",
    "",
    "Not in the old benchmark. The position-only count treats every overlapping peptide in a cluster as an alternative, also when the peptides truly behave alike. Here the alternatives in a stretch are the distinct true patterns of its member peptides, so overlapping peptides with one pattern count once when they share a cluster.",
    "",
    "- **pattern_1k**, **pattern_top**: a pattern is a truth group, i.e. an exact expected profile.",
    "- **pattern_tol_1k**, **pattern_tol_top**: expected profiles within 0.25 log2 (RMS) of one another are one pattern. This is the count for a rule that is meant to let similar overlapping peptides share a cluster.",
    "",
    "Mean log10 over proteins, lower is better.",
    "",
    md_table(summary_overall[
      match(restriction_configs, config), pattern_columns, with = FALSE
    ]), "",
    "Protein by protein, as above:",
    "",
    md_table(potential_paired[, ..pattern_comparison_columns]), "",
    "### How many forms are found, and how many of them are right?",
    "",
    paste0("Counted over the ", forms_summary$proteins[1], " proteins with more than one truth group. ",
           "The restriction changes which clusters end up with members, not how many clusters are asked for."),
    "",
    "- **cl_pairs**: cannot-link pairs the rule gives.",
    "- **present**: true forms, i.e. truth groups of at least two peptides.",
    "- **requested**: clusters asked for (from the dCF labels).",
    "- **found**: forms reported, i.e. clusters with at least two members.",
    "- **recovered**: true forms that one cluster recovers: it holds more than half of the form's peptides and more than half of its members come from the form.",
    "- **recovered_share**: recovered / present. **found_right**: recovered / found.",
    "- **count_error**: mean of |found − present| per protein.",
    "- **profile_R2**: as above, over the same proteins.",
    "",
    md_table(forms_summary), "",
    "Per dataset, as present / found / recovered:",
    "",
    md_table(forms_by_dataset), "",
    "Protein by protein, over all proteins: in how many the first configuration is better / worse than the second, with the sign-test p-value in brackets.",
    "",
    md_table(forms_paired), "",
    "### The switches of the new pipeline",
    "",
    "`new (pipeline)` is the default setting; the three `links` rows change the cannot-link rule, the others use Pearson. **k** is the mean number of clusters used.",
    "",
    md_table(summary_overall[
      grepl("^new", config) & config != "new, unconstrained",
      c("config", "k", "k_error", cluster_columns),
      with = FALSE
    ]), "",
    paste0("Mean number of truth groups with at least two peptides: ",
           formatC(mean(clustering[config == "new (pipeline)", k_true]), format = "f", digits = 2),
           ". `new, true k (reference)` is the new pipeline given that number, to show what a perfect choice of k could gain."),
    "",
    "## 3. Quantification (ProteoMaker)",
    "",
    "Error (log2) of each cluster's quantified profile, as a deviation from its protein's mean profile, against the true deviation of its members. Lower is better.",
    "",
    "The same clustering is quantified four ways, to separate the two changes the new pipeline makes:",
    "",
    "- **all** vs **members**: every peptide contributes, weighted by its membership (old), or only peptides with membership ≥ 0.5 (new).",
    "- **plain** vs **centred**: weighted mean of the intensities (old), or each peptide centred on its own level first (new).",
    "",
    "`quant_all_plain` is the old pipeline's setting and `quant_members_centred` the new one's.",
    "",
    md_table(summary_by_dataset[
      config %in% c("old (pipeline)", "new (pipeline)"),
      .(dataset, clustering = config, quant_all_plain, quant_all_centred,
        quant_members_plain, quant_members_centred)
    ], 3), ""
  )

  # Figure: old (pipeline) -> new (pipeline) per dataset
  metrics_long <- melt(
    pipeline_summary[, .(dataset, config,
             `Assigned\n(higher is better)` = assigned,
             `Pair precision\n(higher is better)` = pair_precision,
             `Pair recall\n(higher is better)` = pair_recall,
             `ARI vs truth groups\n(higher is better)` = ARI,
             `Cannot-links violated\n(lower is better)` = cl_violated,
             `Cluster number error\n(lower is better)` = k_error)],
    id.vars = c("dataset", "config"), variable.name = "metric"
  )
  metrics_long[, dataset := factor(dataset, levels = rev(dataset_levels))]
  metrics_long[, method := ifelse(config == "old (pipeline)", "old", "new")]
  old_and_new <- dcast(metrics_long, dataset + metric ~ method)
  figure <- ggplot() +
    geom_segment(data = old_and_new, aes(x = old, xend = new, y = dataset, yend = dataset),
                 colour = "#c3c2b7", linewidth = 0.7, na.rm = TRUE) +
    geom_point(data = metrics_long, aes(x = value, y = dataset, colour = method),
               size = 3, na.rm = TRUE) +
    scale_colour_manual(values = c(old = old_colour, new = new_colour), breaks = c("old", "new"),
                        labels = c(old = "Old pipeline", new = "New pipeline")) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~metric, nrow = 1, scales = "free_x") +
    labs(title = "Clustering on ProteoMaker, end to end",
         subtitle = "Mean over proteins; each line joins the old and the new pipeline on one dataset") +
    theme_bench() +
    theme(panel.spacing.x = unit(2, "lines"))
  ggsave(file.path(results_dir, "figure_clustering.png"), figure, width = 13.5, height = 3.6, dpi = 200)
  report <- c(report, "![Clustering](figure_clustering.png)", "")
}

# ---------------------------------------------------------------------------
# 4. Real data
# ---------------------------------------------------------------------------
if (!is.null(realdata)) {
  realdata_summary <- realdata[, .(
    complexes = .N,
    peptides = median(as.numeric(peptides)),
    discordant = median(discordant),
    k = median(as.numeric(k)),
    cl_pairs = sum(cl_pairs),
    violated_initial = weighted.mean(violated_initial, cl_pairs, na.rm = TRUE),
    violated_final = weighted.mean(violated_final, cl_pairs, na.rm = TRUE),
    members = median(members),
    nonmember_weight = median(nonmember_weight)
  ), by = pipeline]
  fwrite(realdata_summary, file.path(results_dir, "summary_realdata.csv"))

  per_complex <- dcast(
    realdata, complex ~ pipeline,
    value.var = c("peptides", "discordant", "k", "violated_final")
  )
  per_complex[, complex := substr(complex, 1, 48)]
  setcolorder(per_complex, c("complex", "peptides_old", "peptides_new", "discordant_old",
                             "discordant_new", "k_old", "k_new",
                             "violated_final_old", "violated_final_new"))

  report <- c(
    report,
    "## 4. Real data (H9 complexes)",
    "",
    "No ground truth: this describes what each pipeline does. The one objective check is whether cannot-link pairs end up apart.",
    "",
    "- **peptides**, **discordant**, **k**, **members**, **nonmember_weight**: medians over complexes.",
    "- **violated_initial / violated_final**: share of cannot-link pairs in one cluster before and after the constrained run, pooled over complexes.",
    "- **nonmember_weight**: share of a group's quantification weight that comes from peptides that are not its members.",
    "",
    md_table(realdata_summary), "",
    "### Per complex",
    "",
    md_table(per_complex), ""
  )
}

if (!is.null(restrictions)) {
  has_potential <- "potential_position" %in% names(restrictions)
  has_flexible_cutoff <- "potential_position_1k" %in% names(restrictions)
  if (has_potential) {
    restrictions[, position := log_positive(potential_position)]
    restrictions[,
      position_1k := if (has_flexible_cutoff) {
        log_positive(potential_position_1k)
      } else {
        NA_real_
      }
    ]
    restrictions[, position_top := log_positive(potential_position_top)]
  }
  restriction_summary <- restrictions[, c(
    list(
      complexes = .N,
      requested = sum(k), forms = sum(forms),
      cl_pairs = sum(cl_pairs),
      violated = if (all(is.na(cl_pairs))) NA_real_ else weighted.mean(violated, cl_pairs, na.rm = TRUE),
      members = median(members)
    ),
    if (has_potential) list(position = mean_na(position), position_1k = mean_na(position_1k),
                            position_top = mean_na(position_top))
  ), by = .(pipeline, restriction)]
  # Each restriction against the unconstrained run of its pipeline, complex by complex
  if (has_potential) {
    restrictions_wide <- dcast(
      restrictions, pipeline + complex ~ restriction,
      value.var = c("position", "position_1k", "position_top")
    )
    rules <- setdiff(unique(restrictions$restriction), "none")
    against_unconstrained <- rbindlist(lapply(rules, function(rule) {
      restrictions_wide[, {
        # Sign-test cell: the rule against the unconstrained run
        compare_with_unconstrained <- function(measure) {
          rule_values <- get(paste0(measure, "_", rule))
          unconstrained_values <- get(paste0(measure, "_none"))
          if (!any(is.finite(rule_values) & is.finite(unconstrained_values))) {
            return(NA_character_)
          }
          sign_cell(rule_values, unconstrained_values, lower_is_better = TRUE)
        }
        .(
          restriction = rule,
          position = compare_with_unconstrained("position"),
          position_1k = compare_with_unconstrained("position_1k"),
          position_top = compare_with_unconstrained("position_top")
        )
      }, by = pipeline]
    }))[!is.na(position)]
  }
  report <- c(
    report,
    "### Complexoforms under each restriction",
    "",
    "- **restriction**: `none` is the unconstrained run; the others are cannot-link rules. The old pipeline uses `pearson`; the new pipeline's default is `groups`.",
    "- **requested**: clusters asked for, summed over complexes. **forms**: clusters with at least two members, i.e. the complexoforms reported.",
    "- **cl_pairs**: cannot-link pairs the rule gives. **violated**: share of them in one cluster.",
    "- **members**: share of peptides with a membership ≥ 0.5, median over complexes.",
    "- **position**, **position_1k**, **position_top**: position-only potential proteoforms as in section 2 (mean log10 over complexes, lower is better), counted per protein within each cluster as the old benchmark does for complexes. `position` counts members (membership ≥ 0.5), `position_1k` uses the flexible cutoff (membership ≥ 1/k), `position_top` counts every peptide in its top cluster.",
    "",
    md_table(restriction_summary), ""
  )
  if (has_potential) {
    report <- c(
      report,
      "Complex by complex: in how many a restriction gives fewer / more potential proteoforms than the unconstrained run of the same pipeline, with the sign-test p-value in brackets.",
      "",
      md_table(against_unconstrained), ""
    )
  }
}

writeLines(report, file.path(results_dir, "REPORT.md"))
cat("Wrote", file.path(results_dir, "REPORT.md"), "\n")
