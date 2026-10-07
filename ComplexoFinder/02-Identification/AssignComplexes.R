# Libraries
library(data.table)
library(dynamicTreeCut)
library(limma)

#' Discover complexoforms with peptide-only eBayes prior and direct pool variance
#'
#' Separates the two variance components in the pooled-sibling test:
#'
#' - The **peptide SE** is estimated from each peptide's replicate variance and
#'   shrunk toward a prior built exclusively from peptide-level variances across
#'   all features and conditions, so the prior targets genuinely exchangeable
#'   quantities (replicate noise within individual peptides).
#'
#' - The **pool SE** uses the raw LOO pool variance directly (no eBayes
#'   shrinkage). Because the pool is the leave-one-out mean of all sibling
#'   peptides, its variance naturally captures both measurement noise and
#'   between-peptide profile variance (tau^2), which is exactly the
#'   exchangeability-based H0 variance component required.
#'
#' The Satterthwaite degrees of freedom reflect these two distinct sources:
#' `(d0 + df_i)` for the eBayes-moderated peptide component and `(Nex - 1)`
#' for the direct pool estimate.
#'
#' @inheritParams discover_complexoforms
#' @param cond_adjust_method Character. Multiple-testing correction method used
#'   on the first pass, i.e. across a peptide's per-condition p-values before
#'   the global \code{feature_adjust_method} correction across peptides. One
#'   of \code{stats::p.adjust.methods} (e.g. \code{"BH"} or
#'   \code{"bonferroni"}). Defaults to \code{"BH"}.
#' @param feature_adjust_method Character. Multiple-testing correction method
#'   used on the second pass, i.e. across peptides' first-pass-adjusted
#'   p-values (\code{p_min_bh_feature}), to produce the final \code{p_adj}
#'   column. One of \code{stats::p.adjust.methods} (e.g. \code{"BH"} or
#'   \code{"bonferroni"}). Defaults to \code{"BH"}.
#'
#' @return A list with:
#'
#'   - `peptide_results`: data.table of peptide-level results with discordance
#'     and dCF assignments.
#'   - `condition_results`: data.table of per-condition tests (if requested).
#'
#' @export
discover_complexoforms <- function(
  complex_data,
  intensity_cols,
  group_col = "Gene name",
  peptide_col = "Peptide",
  alpha = 0.05,
  feature_adjust_method = "BH",
  cond_adjust_method = "BH",
  min_total_non_na_frac = 0.40,
  min_conditions = 2,
  min_reps_per_condition = 2,
  deep_split = 2,
  minClusterSize = 1,
  cond_regex = "^(C_[0-9]+)_R_[0-9]+$",
  canonical_label = "dCF0",
  singleton_label = "dCF-1",
  return_condition_tests = FALSE,
  verbose = FALSE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table")
  }
  if (!requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    stop("Please install dynamicTreeCut")
  }
  cond_adjust_method <- match.arg(cond_adjust_method, stats::p.adjust.methods)
  feature_adjust_method <- match.arg(
    feature_adjust_method,
    stats::p.adjust.methods
  )

  DT <- data.table::as.data.table(complex_data)

  # ---- ensure peptide ID ----
  if (!peptide_col %in% names(DT)) {
    needed_cols <- c(
      "Accession",
      "Position in master protein",
      "Modifications in master protein",
      "datatype"
    )
    if (all(needed_cols %in% names(DT))) {
      DT[,
        (peptide_col) := apply(.SD, 1, function(x) paste(x, collapse = "_")),
        .SDcols = needed_cols
      ]
    } else {
      stop(
        "Cannot create peptide identifier column. Provide '",
        peptide_col,
        "' or columns: ",
        paste(needed_cols, collapse = ", ")
      )
    }
  }

  if (!group_col %in% names(DT)) {
    stop("Missing grouping column: ", group_col)
  }
  missing_int <- setdiff(intensity_cols, names(DT))
  if (length(missing_int) > 0) {
    stop("Missing intensity columns: ", paste(missing_int, collapse = ", "))
  }

  # ---- parse conditions ----
  cond_chr <- sub(cond_regex, "\\1", intensity_cols, perl = TRUE)
  if (any(cond_chr == intensity_cols)) {
    alt_regex <- "^[A-Za-z0-9]+_[A-Za-z]+[0-9]+_(D[0-9]+)$"
    cond_chr2 <- sub(alt_regex, "\\1", intensity_cols, perl = TRUE)
    if (any(cond_chr2 == intensity_cols)) {
      stop(
        "Could not parse conditions from intensity_cols.\nTried cond_regex='",
        cond_regex,
        "' and alt_regex='",
        alt_regex,
        "'."
      )
    } else {
      cond_chr <- cond_chr2
      if (verbose) {
        message(
          "Auto-detected H9 format. Using condition pattern: D17, D19, etc."
        )
      }
    }
  }

  cond <- factor(cond_chr)
  cond_levels <- levels(cond)
  K <- length(cond_levels)

  # ---- helpers ----
  .counts_by_condition <- function(Y) {
    counts <- sapply(cond_levels, function(cc) {
      cols_cc <- intensity_cols[cond == cc]
      rowSums(!is.na(Y[, cols_cc, drop = FALSE]))
    })
    as.matrix(counts)
  }

  .double_bh_p <- function(p_row, min_m = 2, method = "BH") {
    p <- p_row[is.finite(p_row)]
    if (length(p) < min_m) {
      return(NA_real_)
    }
    min(stats::p.adjust(p, method = method))
  }

  .cluster_from_deviation <- function(Bd) {
    if (nrow(Bd) < 2) {
      return(setNames(rep(singleton_label, nrow(Bd)), rownames(Bd)))
    }
    Kd <- ncol(Bd)
    n <- nrow(Bd)
    distM <- matrix(0, n, n)
    for (ii in seq_len(n - 1)) {
      xi <- Bd[ii, ]
      for (jj in (ii + 1):n) {
        xj <- Bd[jj, ]
        ok <- !is.na(xi) & !is.na(xj)
        ks <- sum(ok)
        dd <- if (ks >= 2) {
          sqrt(sum((xi[ok] - xj[ok])^2)) * sqrt(Kd / ks)
        } else {
          Inf
        }
        distM[ii, jj] <- dd
        distM[jj, ii] <- dd
      }
    }
    if (any(!is.finite(distM))) {
      fv <- distM[is.finite(distM)]
      distM[!is.finite(distM)] <- (if (length(fv) > 0) max(fv) else 1) * 10
    }
    hc <- stats::hclust(stats::as.dist(distM), method = "ward.D2")
    cl <- if (n >= 3) {
      suppressMessages(dynamicTreeCut::cutreeDynamic(
        dendro = hc,
        distM = distM,
        minClusterSize = minClusterSize,
        method = "hybrid",
        deepSplit = deep_split,
        verbose = 0
      ))
    } else {
      seq_len(n)
    }
    cl <- as.integer(cl)
    names(cl) <- rownames(Bd)
    cl_dt <- data.table::data.table(peptide = names(cl), cluster_id = cl)
    sizes <- cl_dt[, .N, by = cluster_id]
    good_ids <- sort(unique(sizes[cluster_id > 0 & N >= 2, cluster_id]))
    cl_dt[, label := singleton_label]
    if (length(good_ids) > 0) {
      map <- data.table::data.table(
        cluster_id = good_ids,
        label = paste0("dCF", seq_along(good_ids))
      )
      cl_dt <- merge(
        cl_dt,
        map,
        by = "cluster_id",
        all.x = TRUE,
        sort = FALSE,
        suffixes = c("", ".m")
      )
      cl_dt[!is.na(label.m), label := label.m]
      cl_dt[, label.m := NULL]
    }
    stats::setNames(cl_dt$label, cl_dt$peptide)
  }

  if (verbose) {
    message(
      "Running PEPTIDE-PRIOR eBayes pooled-sibling analysis on ",
      nrow(DT),
      " peptides from ",
      data.table::uniqueN(DT[[group_col]]),
      " proteins"
    )
  }

  outg <- DT[, c(group_col, peptide_col, intensity_cols), with = FALSE]
  outg[, discordant := FALSE]
  outg[, dCF := canonical_label]

  Y <- as.matrix(DT[, ..intensity_cols])
  storage.mode(Y) <- "double"

  # ---- filter peptides ----
  n_samples <- ncol(Y)
  min_total_non_na <- ceiling(n_samples * min_total_non_na_frac)
  counts <- .counts_by_condition(Y)
  keep <- (rowSums(!is.na(Y)) >= min_total_non_na) &
    (rowSums(counts >= min_reps_per_condition) >= min_conditions)

  if (verbose) {
    message(
      "Filtering: n_samples=",
      n_samples,
      ", min_total_non_na=",
      min_total_non_na,
      ", K (conditions)=",
      K
    )
    message(
      "  Peptides passing total filter: ",
      sum(rowSums(!is.na(Y)) >= min_total_non_na)
    )
    message(
      "  Peptides passing condition filter: ",
      sum(rowSums(counts >= min_reps_per_condition) >= min_conditions)
    )
    message("  Peptides passing both filters: ", sum(keep))
  }

  if (sum(keep) < 2) {
    if (verbose) {
      message(
        "Insufficient peptides passed filtering. Returning all as canonical."
      )
    }
    outg[, n_cond_tested := NA_integer_]
    outg[, p_min_bh_feature := NA_real_]
    outg[, p_adj := NA_real_]
    return(list(
      peptide_results = outg[,
        c(
          group_col,
          peptide_col,
          intensity_cols,
          "n_cond_tested",
          "p_min_bh_feature",
          "p_adj",
          "discordant",
          "dCF"
        ),
        with = FALSE
      ],
      condition_results = data.table::data.table()
    ))
  }

  Yk <- Y[keep, , drop = FALSE]
  pepK <- DT[[peptide_col]][keep]
  P <- nrow(Yk)

  mk <- function(fill = NA_real_) {
    matrix(fill, nrow = P, ncol = K, dimnames = list(pepK, cond_levels))
  }
  mean_i <- mk()
  var_i <- mk()
  n_i <- mk(0L)
  pool_mean <- mk()
  pool_var <- mk()
  pool_n <- mk(0L)
  t_mat <- mk()
  df_mat <- mk()
  p_mat <- mk()
  fc_mat <- mk()

  # ---- Pass 1: per-condition LOO statistics ----
  for (kk in seq_along(cond_levels)) {
    cc <- cond_levels[kk]
    cols_cc <- intensity_cols[cond == cc]
    X <- Yk[, cols_cc, drop = FALSE]

    ni <- rowSums(!is.na(X))
    si <- rowSums(X, na.rm = TRUE)
    ssi <- rowSums(X * X, na.rm = TRUE)
    mi <- si / ni
    vi <- (ssi - (si * si) / ni) / pmax(ni - 1, 1)

    Ntot <- sum(ni)
    Stot <- sum(si, na.rm = TRUE)
    SStot <- sum(ssi, na.rm = TRUE)
    Nex <- Ntot - ni
    Sex <- Stot - si
    SSex <- SStot - ssi
    mex <- Sex / Nex
    vex <- (SSex - (Sex * Sex) / Nex) / pmax(Nex - 1, 1)

    n_i[, kk] <- as.integer(ni)
    mean_i[, kk] <- mi
    var_i[, kk] <- vi
    pool_n[, kk] <- as.integer(Nex)
    pool_mean[, kk] <- mex
    pool_var[, kk] <- vex
    fc_mat[, kk] <- mi - mex
  }

  # ---- EMPIRICAL BAYES VARIANCE MODERATION (peptide variances only) ----
  # Prior is estimated exclusively from peptide-level replicate variances,
  # which are the only exchangeable quantities in this context.
  # Pool variances (sigma^2 + tau^2) are a different kind of quantity and
  # are excluded to avoid contaminating the peptide prior.
  pep_vars <- as.vector(var_i)
  pep_df <- as.vector(pmax(n_i - 1, 1))
  ok_pep <- is.finite(pep_df) & pep_df > 0 & is.finite(pep_vars) & pep_vars > 0

  eb_method <- "moM"
  d0 <- 4
  s0_sq <- stats::median(pep_vars[ok_pep], na.rm = TRUE)

  if (sum(ok_pep) >= 10) {
    if (requireNamespace("limma", quietly = TRUE)) {
      eb_fit <- tryCatch(
        limma::fitFDist(pep_vars[ok_pep], df1 = pep_df[ok_pep]),
        error = function(e) NULL
      )
      if (!is.null(eb_fit)) {
        d0 <- eb_fit$df1
        s0_sq <- eb_fit$scale
        eb_method <- "limma::fitFDist (peptide-only)"
      } else {
        mv <- mean(pep_vars[ok_pep])
        vv <- stats::var(pep_vars[ok_pep])
        if (vv > 0 && mv > 0) {
          d0 <- max(2 * mv^2 / vv + 4, 1)
          s0_sq <- mv * (d0 - 2) / d0
        }
      }
    } else {
      mv <- mean(pep_vars[ok_pep])
      vv <- stats::var(pep_vars[ok_pep])
      if (vv > 0 && mv > 0) {
        d0 <- max(2 * mv^2 / vv + 4, 1)
        s0_sq <- mv * (d0 - 2) / d0
      }
    }
  }

  if (!is.numeric(d0) || length(d0) != 1 || !is.finite(d0) || d0 <= 0) {
    mv <- mean(pep_vars[ok_pep], na.rm = TRUE)
    vv <- stats::var(pep_vars[ok_pep], na.rm = TRUE)
    d0 <- if (vv > 0 && mv > 0) max(2 * mv^2 / vv + 4, 1) else 4
    s0_sq <- if (vv > 0 && mv > 0) mv * (d0 - 2) / d0 else mv
    eb_method <- "moM"
  }
  if (!is.numeric(s0_sq) || !is.finite(s0_sq) || s0_sq <= 0) {
    s0_sq <- 1e-6
  }

  if (verbose) {
    message("Empirical Bayes variance moderation (peptide prior only):")
    message("  Method    = ", eb_method)
    message("  Prior df  (d0)    = ", round(d0, 2))
    message("  Prior var (s0^2)  = ", round(s0_sq, 4))
  }

  # ---- Pass 2: moderated test statistics ----
  # Peptide variance is eBayes-shrunk toward the peptide-only prior.
  # Pool variance is used as a direct sample estimate (no shrinkage); its
  # df = Nex - 1 enters the Satterthwaite formula directly.
  for (kk in seq_along(cond_levels)) {
    ok <- is.finite(mean_i[, kk]) &
      is.finite(pool_mean[, kk]) &
      (n_i[, kk] >= min_reps_per_condition) &
      (pool_n[, kk] >= 2) &
      is.finite(var_i[, kk]) &
      var_i[, kk] >= 0 &
      is.finite(pool_var[, kk]) &
      pool_var[, kk] >= 0

    if (any(ok)) {
      df_i <- pmax(n_i[ok, kk] - 1, 1)

      # Moderate peptide variance toward the peptide-only prior
      var_i_mod <- (d0 * s0_sq + df_i * var_i[ok, kk]) / (d0 + df_i)

      # Pool SE: raw estimate — no eBayes shrinkage.
      # var_pool captures sigma^2 + tau^2 under H0 (exchangeability), giving
      # the correct denominator variance for the pool mean.
      df_pool <- pmax(pool_n[ok, kk] - 1, 1)
      a <- var_i_mod / n_i[ok, kk]
      b <- pool_var[ok, kk] / pool_n[ok, kk]

      se2 <- a + b
      tval <- (mean_i[ok, kk] - pool_mean[ok, kk]) / sqrt(se2)

      # Satterthwaite df: eBayes effective df for peptide, direct df for pool
      df_satt <- (a + b)^2 / ((a^2) / (d0 + df_i) + (b^2) / df_pool)

      t_mat[ok, kk] <- tval
      df_mat[ok, kk] <- df_satt
      p_mat[ok, kk] <- 2 * stats::pt(-abs(tval), df = df_satt)
    }
  }

  # ---- feature-level BH over conditions, then global BH ----
  n_cond_tested <- apply(p_mat, 1, function(x) sum(is.finite(x)))
  p_min_bh_feature <- apply(
    p_mat,
    1,
    .double_bh_p,
    min_m = min_conditions,
    method = cond_adjust_method
  )

  p_adj <- rep(NA_real_, length(p_min_bh_feature))
  ok_s <- is.finite(p_min_bh_feature)
  if (any(ok_s)) {
    p_adj[ok_s] <- stats::p.adjust(
      p_min_bh_feature[ok_s],
      method = feature_adjust_method
    )
  }

  disc <- ok_s & is.finite(p_adj) & (p_adj < alpha)

  idx <- which(keep)
  outg[, n_cond_tested := NA_integer_]
  outg[, p_min_bh_feature := NA_real_]
  outg[, p_adj := NA_real_]
  data.table::set(
    outg,
    i = idx,
    j = "n_cond_tested",
    value = as.integer(n_cond_tested)
  )
  data.table::set(
    outg,
    i = idx,
    j = "p_min_bh_feature",
    value = p_min_bh_feature
  )
  data.table::set(outg, i = idx, j = "p_adj", value = p_adj)
  outg[idx, discordant := disc]

  if (verbose) {
    message(sum(disc, na.rm = TRUE), " peptides marked as discordant")
  }

  # ---- clustering on deviation profiles ----
  deviation <- mean_i - pool_mean
  rownames(deviation) <- pepK
  disc_peps <- pepK[disc]
  outg[, dCF := canonical_label]

  if (length(disc_peps) == 1) {
    outg[get(peptide_col) %in% disc_peps, dCF := singleton_label]
  } else if (length(disc_peps) > 1) {
    Bd <- deviation[disc_peps, , drop = FALSE]
    labs <- .cluster_from_deviation(Bd)
    outg[get(peptide_col) %in% names(labs), dCF := labs[get(peptide_col)]]
  }

  # ---- condition-level table ----
  cond_res <- data.table::data.table()
  if (return_condition_tests) {
    cond_res <- data.table::as.data.table(as.table(p_mat))
    data.table::setnames(cond_res, c("peptide", "condition", "p_value"))
    cond_res[, t_stat := as.vector(t_mat)]
    cond_res[, df := as.vector(df_mat)]
    cond_res[, deviation := as.vector(deviation)]
    cond_res[, fold_change := as.vector(fc_mat)]
    cond_res[, mean_peptide := as.vector(mean_i)]
    cond_res[, mean_pool := as.vector(pool_mean)]
    cond_res[, n_peptide := as.vector(n_i)]
    cond_res[, n_pool := as.vector(pool_n)]
    cond_res <- cond_res[is.finite(p_value)]
  }

  list(
    peptide_results = outg[,
      c(
        group_col,
        peptide_col,
        intensity_cols,
        "n_cond_tested",
        "p_min_bh_feature",
        "p_adj",
        "discordant",
        "dCF"
      ),
      with = FALSE
    ],
    condition_results = cond_res
  )
}
