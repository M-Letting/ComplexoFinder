# ============================================================================
# 01-Discovery.R
# ============================================================================
# Discordant-peptide discovery, old against new, on the two benchmarks with a
# ground truth (see BenchmarkData.R):
#
#   old  discover_complexoforms() of ComplexoFinder, once per protein, with
#        the pipeline's settings
#   new  discover_complexoforms() of NewImplementation, one call per dataset
#   new_global_bh  the same with correction = "global_bh": Benjamini-Hochberg
#        across all peptides only, without the Bonferroni step within each
#        protein
#
# For SWATH-MS the saved results of the published tools (PeCorA, ProteoForge,
# COPF; Benchmark/01-FindDCF/data/results, not re-run here) are scored the
# same way, as a reference.
#
# Writes results/discovery_metrics.csv and results/discovery_calls.rds.
#
# Run from the project root (~10 min):
#   Rscript NewImplementation/Benchmark/01-Discovery.R
# ============================================================================

suppressMessages(library(data.table))
setwd(here::here())
source("NewImplementation/Benchmark/BenchmarkData.R")
source("NewImplementation/Benchmark/BenchmarkFunctions.R")

results_dir <- "NewImplementation/Benchmark/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

old_pipeline <- load_old_implementation()
new_pipeline <- load_new_implementation()

benchmarks <- list(
  "SWATH-MS" = list(datasets = swath_datasets, loader = load_swath_dataset),
  "ProteoMaker" = list(datasets = proteomaker_datasets, loader = load_proteomaker_dataset)
)

published_results_dir <- "Benchmark/01-FindDCF/data/results"

#' Read the saved calls of a published tool on a SWATH-MS dataset
#'
#' The score is the tool's own p-value; a peptide is called when its adjusted
#' p-value is below 0.05. COPF scores proteins, not peptides, so it has no
#' peptide-level calls.
#'
#' @param tool "PeCorA", "ProteoForge" or "COPF".
#' @param dataset_name One of `swath_datasets`.
#'
#' @return data.table with `peptide_id`, `score` and `called`, or `NULL` if
#'   there is no saved result.
read_published <- function(tool, dataset_name) {
  result_file <- file.path(
    published_results_dir,
    paste0(tool, "_", dataset_name, "_result.feather")
  )
  if (!file.exists(result_file)) return(NULL)
  tool_results <- as.data.table(as.data.frame(arrow::read_feather(result_file)))
  switch(
    tool,
    PeCorA = unique(tool_results[, .(
      peptide_id, score = pvalue, called = adj_pval < 0.05
    )]),
    ProteoForge = unique(tool_results[, .(
      peptide_id, score = pval, called = adj_pval < 0.05
    )]),
    COPF = unique(tool_results[, .(
      peptide_id = id, score = proteoform_score_pval, called = NA
    )])
  )
}

metrics <- list()
calls <- list()

for (benchmark_name in names(benchmarks)) {
  benchmark <- benchmarks[[benchmark_name]]
  for (dataset_name in benchmark$datasets) {
    dataset <- benchmark$loader(dataset_name)
    cat(sprintf(
      "\n[%s / %s] %d peptides, %d assemblies\n",
      benchmark_name, dataset_name,
      nrow(dataset$data), uniqueN(dataset$data[[dataset$assembly_column]])
    ))

    calls_by_method <- list()
    start_time <- Sys.time()
    calls_by_method$old <- run_discovery_old(old_pipeline, dataset)
    old_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    start_time <- Sys.time()
    calls_by_method$new <- run_discovery_new(new_pipeline, dataset)$calls
    new_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    calls_by_method$new_global_bh <- run_discovery_new(
      new_pipeline, dataset, correction = "global_bh"
    )$calls
    seconds <- c(old = old_seconds, new = new_seconds, new_global_bh = new_seconds)

    if (benchmark_name == "SWATH-MS") {
      for (tool in c("PeCorA", "ProteoForge", "COPF")) {
        tool_calls <- read_published(tool, dataset_name)
        if (is.null(tool_calls)) next
        tool_calls[, peptide_id := as.character(peptide_id)]
        # Peptide ids are unique across proteins in this benchmark
        tool_calls[,
          group_id := dataset$truth$group_id[
            match(peptide_id, as.character(dataset$truth$peptide_id))
          ]
        ]
        calls_by_method[[tool]] <- tool_calls
        seconds[tool] <- NA
      }
    }

    for (method in names(calls_by_method)) {
      method_calls <- calls_by_method[[method]]
      method_metrics <- evaluate_discovery(method_calls, dataset$truth)
      # A tool without peptide-level calls has no false-alarm rate or power
      if (all(is.na(method_calls$called))) {
        method_metrics[, c("fpr", "tpr", "fdp") := NA_real_]
      }
      metrics[[length(metrics) + 1]] <- cbind(
        data.table(benchmark = benchmark_name, dataset = dataset_name, method = method),
        method_metrics,
        seconds = seconds[[method]]
      )
      calls[[length(calls) + 1]] <- cbind(
        data.table(benchmark = benchmark_name, dataset = dataset_name, method = method),
        merge(
          dataset$truth[, .(group_id, peptide_id, truth_discordant)],
          method_calls[, .(group_id, peptide_id, score, called)],
          by = c("group_id", "peptide_id"),
          all.x = TRUE
        )
      )
      cat(sprintf(
        "  %-14s coverage %.2f  AUROC %.3f  within %.3f  false alarms %.3f  power %.3f\n",
        method, method_metrics$coverage, method_metrics$auroc, method_metrics$within,
        method_metrics$fpr, method_metrics$tpr
      ))
    }
  }
}

fwrite(rbindlist(metrics), file.path(results_dir, "discovery_metrics.csv"))
saveRDS(rbindlist(calls, fill = TRUE), file.path(results_dir, "discovery_calls.rds"))
cat("\nWrote", file.path(results_dir, "discovery_metrics.csv"), "\n")
