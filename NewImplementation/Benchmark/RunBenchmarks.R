# ============================================================================
# RunBenchmarks.R
# ============================================================================
# Runs the whole benchmark suite: the old implementation (ComplexoFinder/)
# against the new one (NewImplementation/), then the report.
#
#   01-Discovery.R    discordant peptides, SWATH-MS and ProteoMaker   ~2 min
#   02-Clustering.R   clustering and quantification, ProteoMaker     ~35 min
#   02b-CannotLinks.R  the cannot-link rules pair by pair, ProteoMaker  ~3 min
#   03-RealData.R     the H9 complexes (needs the new pipeline's
#                     results in NewImplementation/10-Results/H9)     ~15 min
#   04-Report.R       results/REPORT.md, summary tables and figures
#
# Each script runs in its own R session and can also be run on its own; the
# report uses whichever results exist.
#
# Run from the project root:
#   Rscript NewImplementation/Benchmark/RunBenchmarks.R
# ============================================================================

setwd(here::here())

scripts <- c("01-Discovery.R", "02-Clustering.R", "02b-CannotLinks.R", "03-RealData.R", "04-Report.R")
rscript_command <- file.path(R.home("bin"), "Rscript")

for (script in scripts) {
  cat(paste(rep("#", 80), collapse = ""), "\n")
  cat("Running", script, "\n")
  start_time <- Sys.time()
  exit_status <- system2(rscript_command, file.path("NewImplementation/Benchmark", script))
  cat(sprintf(
    "%s finished in %.1f min%s\n",
    script,
    as.numeric(difftime(Sys.time(), start_time, units = "mins")),
    if (exit_status != 0) " WITH ERRORS" else ""
  ))
}
