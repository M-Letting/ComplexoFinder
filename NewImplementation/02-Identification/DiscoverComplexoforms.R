# Libraries
library(data.table)
library(limma)

# ============================================================================
# DISCORDANT-PEPTIDE DISCOVERY (DCF_ortho_relaxed)
# ============================================================================
# A peptide is discordant when its PROFILE across conditions differs from the
# profile shared by the largest self-consistent group of its assembly's
# peptides (the coherent core). The method is described in
# Misc/DCF_ortho_relaxed_Methodology.md, which writes se^2 for the squared
# standard error of a condition mean and tau^2 for the between-peptide
# variance.
#
# The entry point, discover_complexoforms(), is called ONCE on the whole
# multi-assembly table (every complex stacked, one column naming the complex):
# the variance priors, the between-peptide variance and the final
# Benjamini-Hochberg step are estimated across all assemblies.
#
# Nothing is imputed: peptides are kept or dropped, never filled.
# ============================================================================

#' Sort condition labels by the number they contain
#'
#' @param conditions Character vector of condition labels.
#' @return The labels sorted by their digits (\code{"D9"} before
#'   \code{"D10"}), or alphabetically if a label has no digits.
sort_conditions_nat <- function(conditions) {
  condition_numbers <- suppressWarnings(as.numeric(gsub("[^0-9]", "", conditions)))
  if (any(is.na(condition_numbers))) {
    conditions[order(conditions)]
  } else {
    conditions[order(condition_numbers)]
  }
}

#' Parse conditions and apply the completeness filter
#'
#' A peptide is kept when at least \code{min_fraction_observed} of its values
#' are present and at least \code{min_conditions} conditions hold
#' \code{min_replicates_per_condition} values. An assembly needs two kept
#' peptides, so that there is another peptide to compare with.
#'
#' @param data data.table (or coercible) with one row per assembly and peptide.
#' @param intensity_columns Character vector of intensity column names.
#' @param assembly_column Character. Column identifying the assembly.
#' @param peptide_column Character. Column with peptide identifiers.
#' @param condition_regex Regex with one capture group extracting the condition
#'   label from an intensity column name.
#' @param min_fraction_observed Numeric in [0, 1]. Minimum fraction of
#'   non-missing intensities.
#' @param min_conditions Integer. Minimum number of conditions that must meet
#'   \code{min_replicates_per_condition}.
#' @param min_replicates_per_condition Integer. Minimum non-missing replicates
#'   for a condition to count.
#' @param datatype_column Character or \code{NULL}. Column with the datatype.
#' @param verbose Logical.
#'
#' @return List with \code{tested} (logical, one per row of \code{data}: the
#'   peptide passes the filter), \code{intensities} (matrix), \code{assembly},
#'   \code{peptide} and \code{datatype} (the last four for the tested peptides
#'   only), \code{column_conditions} (condition label of each intensity column)
#'   and \code{condition_levels} (the conditions in order).
prepare_peptides <- function(
  data,
  intensity_columns,
  assembly_column,
  peptide_column,
  condition_regex,
  min_fraction_observed = 0,
  min_conditions = 2,
  min_replicates_per_condition = 1,
  datatype_column = NULL,
  verbose = FALSE
) {
  peptide_table <- as.data.table(data)
  for (column in c(assembly_column, peptide_column, datatype_column)) {
    if (!column %in% names(peptide_table)) stop("Missing column: ", column)
  }
  missing_intensity_columns <- setdiff(intensity_columns, names(peptide_table))
  if (length(missing_intensity_columns)) {
    stop(
      "Missing intensity columns: ",
      paste(missing_intensity_columns, collapse = ", ")
    )
  }

  column_conditions <- sub(condition_regex, "\\1", intensity_columns, perl = TRUE)
  if (any(column_conditions == intensity_columns)) {
    stop(
      "Could not parse conditions from intensity_columns with condition_regex='",
      condition_regex,
      "'"
    )
  }
  condition_levels <- sort_conditions_nat(unique(column_conditions))

  all_intensities <- as.matrix(peptide_table[, intensity_columns, with = FALSE])
  storage.mode(all_intensities) <- "double"

  # Observed replicates per peptide and condition
  replicates_observed <- vapply(
    condition_levels,
    function(condition) {
      rowSums(!is.na(all_intensities[, column_conditions == condition, drop = FALSE]))
    },
    numeric(nrow(all_intensities))
  )
  replicates_observed <- matrix(replicates_observed, nrow = nrow(all_intensities))

  tested <- (rowSums(!is.na(all_intensities)) >=
    ceiling(ncol(all_intensities) * min_fraction_observed)) &
    (rowSums(replicates_observed >= min_replicates_per_condition) >= min_conditions)

  # An assembly needs at least two tested peptides
  assembly <- as.character(peptide_table[[assembly_column]])
  peptides_per_assembly <- table(assembly[tested])
  tested <- tested &
    (assembly %in% names(peptides_per_assembly)[peptides_per_assembly >= 2])

  if (verbose) {
    message(
      "Kept ",
      sum(tested),
      " / ",
      nrow(peptide_table),
      " peptides in ",
      length(unique(assembly[tested])),
      " assemblies"
    )
  }

  list(
    tested = tested,
    intensities = all_intensities[tested, , drop = FALSE],
    assembly = assembly[tested],
    peptide = as.character(peptide_table[[peptide_column]])[tested],
    datatype = if (!is.null(datatype_column)) {
      as.character(peptide_table[[datatype_column]])[tested]
    } else {
      NULL
    },
    column_conditions = column_conditions,
    condition_levels = condition_levels
  )
}

#' Mean, variance and replicate count per peptide and condition
#'
#' @param intensities Numeric matrix, peptides x samples.
#' @param column_conditions Character vector, condition label of each column.
#' @param condition_levels Character vector of the conditions, in order.
#' @param min_replicates_per_condition Integer. A condition mean is \code{NA}
#'   when the peptide has fewer observed replicates in the condition.
#'
#' @return List of three matrices, peptides x conditions: \code{mean},
#'   \code{variance} (\code{NA} with fewer than two replicates) and
#'   \code{replicates} (number of observed replicates).
condition_stats <- function(
  intensities,
  column_conditions,
  condition_levels,
  min_replicates_per_condition = 1
) {
  n_peptides <- nrow(intensities)
  n_conditions <- length(condition_levels)
  filled_matrix <- function(value = NA_real_) {
    matrix(value, n_peptides, n_conditions, dimnames = list(NULL, condition_levels))
  }
  condition_mean <- filled_matrix()
  condition_variance <- filled_matrix()
  replicates <- filled_matrix(0)
  for (condition_index in seq_len(n_conditions)) {
    condition_intensities <- intensities[,
      column_conditions == condition_levels[condition_index],
      drop = FALSE
    ]
    n_observed <- rowSums(!is.na(condition_intensities))
    value_sum <- rowSums(condition_intensities, na.rm = TRUE)
    squared_sum <- rowSums(condition_intensities * condition_intensities, na.rm = TRUE)
    replicates[, condition_index] <- n_observed
    condition_mean[, condition_index] <- ifelse(
      n_observed > 0, value_sum / n_observed, NA_real_
    )
    condition_variance[, condition_index] <- ifelse(
      n_observed > 1,
      (squared_sum - (value_sum * value_sum) / pmax(n_observed, 1)) /
        pmax(n_observed - 1, 1),
      NA_real_
    )
  }
  condition_mean[replicates < min_replicates_per_condition] <- NA_real_
  list(mean = condition_mean, variance = condition_variance, replicates = replicates)
}

#' Pool the replicate variance across conditions and moderate it
#'
#' A peptide's replicate variance is the average of its condition variances,
#' weighted by their degrees of freedom. It is then shrunk towards a prior
#' fitted over the peptides of its stratum (empirical Bayes,
#' \code{limma::fitFDist()}; with too few peptides, a prior at the median
#' variance with 4 degrees of freedom). A peptide with no condition holding two
#' replicates takes the prior variance.
#'
#' @param condition_variance Numeric matrix, peptides x conditions.
#' @param replicates Numeric matrix, peptides x conditions: observed
#'   replicates.
#' @param prior_stratum Optional character vector, one per peptide (e.g. the
#'   datatype): one prior is fitted per stratum instead of one for all.
#' @param min_peptides_for_prior Integer. Number of peptides a stratum needs
#'   for its prior to be fitted with \code{limma::fitFDist()}.
#' @param verbose Logical. Print the fitted priors.
#'
#' @return List with \code{pooled_variance} (before moderation),
#'   \code{moderated_variance}, \code{degrees_of_freedom} (of the moderated
#'   variance) and \code{prior_table} (one row per stratum: \code{stratum},
#'   prior degrees of freedom \code{d0}, prior variance \code{s0_sq},
#'   \code{eb_method} and number of peptides \code{n}).
moderate_variance <- function(
  condition_variance,
  replicates,
  prior_stratum = NULL,
  min_peptides_for_prior = 10,
  verbose = FALSE
) {
  condition_degrees_of_freedom <- pmax(replicates - 1, 0)
  usable <- is.finite(condition_variance) & condition_degrees_of_freedom > 0
  weighted_variance_sum <- rowSums(
    ifelse(usable, condition_degrees_of_freedom * condition_variance, 0),
    na.rm = TRUE
  )
  pooled_degrees_of_freedom <- rowSums(
    ifelse(usable, condition_degrees_of_freedom, 0),
    na.rm = TRUE
  )
  pooled_variance <- ifelse(
    pooled_degrees_of_freedom > 0,
    weighted_variance_sum / pooled_degrees_of_freedom,
    NA_real_
  )
  pooled_degrees_of_freedom <- pmax(pooled_degrees_of_freedom, 0)

  stratum <- if (is.null(prior_stratum)) {
    rep("GLOBAL", length(pooled_variance))
  } else {
    prior_stratum
  }
  moderated_variance <- rep(NA_real_, length(pooled_variance))
  moderated_degrees_of_freedom <- rep(NA_real_, length(pooled_variance))
  prior_rows <- list()

  for (stratum_name in unique(stratum)) {
    in_stratum <- which(stratum == stratum_name)
    informative <- in_stratum[
      is.finite(pooled_variance[in_stratum]) &
        pooled_variance[in_stratum] > 0 &
        pooled_degrees_of_freedom[in_stratum] > 0
    ]
    prior_degrees_of_freedom <- 4
    prior_variance <- if (length(informative)) {
      stats::median(pooled_variance[informative])
    } else {
      1e-6
    }
    prior_method <- "moM"
    if (length(informative) >= min_peptides_for_prior) {
      prior_fit <- tryCatch(
        limma::fitFDist(
          pooled_variance[informative],
          df1 = pooled_degrees_of_freedom[informative]
        ),
        error = function(e) NULL
      )
      if (!is.null(prior_fit) && is.finite(prior_fit$df2) && prior_fit$df2 > 0) {
        prior_degrees_of_freedom <- min(prior_fit$df2, 1e6)
        prior_variance <- prior_fit$scale
        prior_method <- "limma::fitFDist"
      }
    }
    if (!is.finite(prior_variance) || prior_variance <= 0) prior_variance <- 1e-6
    moderated_variance[in_stratum] <- (prior_degrees_of_freedom * prior_variance +
      pooled_degrees_of_freedom[in_stratum] *
        ifelse(
          is.finite(pooled_variance[in_stratum]),
          pooled_variance[in_stratum],
          prior_variance
        )) /
      (prior_degrees_of_freedom + pooled_degrees_of_freedom[in_stratum])
    moderated_degrees_of_freedom[in_stratum] <- prior_degrees_of_freedom +
      pooled_degrees_of_freedom[in_stratum]
    prior_rows[[stratum_name]] <- data.table(
      stratum = stratum_name,
      d0 = prior_degrees_of_freedom,
      s0_sq = prior_variance,
      eb_method = prior_method,
      n = length(informative)
    )
  }
  if (verbose) print(rbindlist(prior_rows))
  list(
    pooled_variance = pooled_variance,
    moderated_variance = moderated_variance,
    degrees_of_freedom = moderated_degrees_of_freedom,
    prior_table = rbindlist(prior_rows)
  )
}

#' Split each peptide into its offset and its centered profile
#'
#' The offset is the peptide's mean over the conditions it is observed in (its
#' own abundance level). The centered profile is what remains: how the peptide
#' changes across conditions.
#'
#' @param condition_mean Numeric matrix, peptides x conditions.
#'
#' @return List with \code{offset} (per peptide), \code{profile} (matrix, the
#'   condition means minus the offset) and \code{conditions_observed} (number
#'   of conditions with a value, per peptide).
split_offset_profile <- function(condition_mean) {
  conditions_observed <- rowSums(is.finite(condition_mean))
  peptide_mean <- rowMeans(condition_mean, na.rm = TRUE)
  peptide_mean[conditions_observed == 0] <- NA_real_
  list(
    offset = peptide_mean,
    profile = condition_mean - peptide_mean,
    conditions_observed = conditions_observed
  )
}

#' Standardized pairwise profile distance within one assembly
#'
#' The distance between two peptides is the root mean square, over the
#' conditions both are observed in, of the difference between their centered
#' profiles divided by its standard deviation (the square root of the two
#' peptides' variances added). Two peptides with the same profile are about
#' one unit apart.
#'
#' @param centered_profiles Numeric matrix, peptides x conditions.
#' @param profile_variance Numeric matrix, peptides x conditions: variance of
#'   each profile value.
#' @param min_shared_conditions Integer. Pairs with fewer jointly observed
#'   conditions get \code{NA}.
#'
#' @return Numeric matrix, peptides x peptides.
.pairwise_profile_distance <- function(
  centered_profiles,
  profile_variance,
  min_shared_conditions = 2
) {
  n_peptides <- nrow(centered_profiles)
  standardized_square_sum <- matrix(0, n_peptides, n_peptides)
  shared_conditions <- matrix(0L, n_peptides, n_peptides)
  for (condition_index in seq_len(ncol(centered_profiles))) {
    profile_value <- centered_profiles[, condition_index]
    variance <- profile_variance[, condition_index]
    observed <- is.finite(profile_value) & is.finite(variance) & variance > 0
    if (sum(observed) < 2) next
    difference <- outer(profile_value, profile_value, "-")
    difference_variance <- outer(variance, variance, "+")
    comparable <- outer(observed, observed, "&") &
      is.finite(difference_variance) &
      difference_variance > 0
    standardized_square_sum <- standardized_square_sum +
      ifelse(comparable, difference * difference / difference_variance, 0)
    shared_conditions <- shared_conditions + comparable
  }
  distance <- ifelse(
    shared_conditions >= min_shared_conditions,
    sqrt(standardized_square_sum / pmax(shared_conditions, 1)),
    NA_real_
  )
  diag(distance) <- 0
  distance
}

#' Hierarchical partition of one assembly's peptides by profile
#'
#' Clusters the peptides hierarchically on the standardized profile distance
#' (\code{.pairwise_profile_distance()}) and cuts the tree. Pairs whose
#' distance cannot be measured (too few shared conditions) are treated as far
#' apart, so they cannot pull clusters together.
#'
#' The tree is cut in one of two ways:
#' \describe{
#'   \item{\code{dynamic = FALSE}}{Average linkage, cut at the fixed height
#'     \code{cut_height}.}
#'   \item{\code{dynamic = TRUE}}{Ward.D2 linkage, cut by
#'     \code{dynamicTreeCut::cutreeDynamic()} with a minimum cluster size of
#'     one, so a peptide is never merged into a cluster to reach a size. With
#'     fewer than three peptides there is no tree to cut dynamically, and the
#'     fixed height is used.}
#' }
#'
#' @param centered_profiles Numeric matrix, peptides x conditions.
#' @param profile_variance Numeric matrix, peptides x conditions: variance of
#'   each profile value.
#' @param cut_height Numeric. Height of the fixed cut.
#' @param min_shared_conditions Integer. Conditions two peptides must share for
#'   their distance to be measured.
#' @param dynamic Logical. Use the dynamic tree cut.
#' @param deep_split Integer 0-4. Sensitivity of the dynamic tree cut to
#'   splitting clusters (\code{deepSplit}); higher gives more, smaller
#'   clusters.
#'
#' @return Integer vector of cluster numbers, one per peptide; \code{NA} for a
#'   peptide with no measurable distance to any other.
.partition_profiles <- function(
  centered_profiles,
  profile_variance,
  cut_height,
  min_shared_conditions = 2,
  dynamic = FALSE,
  deep_split = 2
) {
  cluster <- rep(NA_integer_, nrow(centered_profiles))
  distance <- .pairwise_profile_distance(
    centered_profiles,
    profile_variance,
    min_shared_conditions
  )
  measurable <- rowSums(is.finite(distance)) > 1
  if (sum(measurable) < 2) {
    return(cluster)
  }
  measurable_distance <- distance[measurable, measurable, drop = FALSE]
  measurable_distance[!is.finite(measurable_distance)] <- max(
    c(measurable_distance[is.finite(measurable_distance)], cut_height),
    na.rm = TRUE
  ) * 2

  dynamic <- dynamic && sum(measurable) >= 3
  tree <- tryCatch(
    stats::hclust(
      stats::as.dist(measurable_distance),
      method = if (dynamic) "ward.D2" else "average"
    ),
    error = function(e) NULL
  )
  if (is.null(tree)) {
    return(cluster)
  }
  if (dynamic) {
    tree_cluster <- as.integer(suppressMessages(dynamicTreeCut::cutreeDynamic(
      dendro = tree,
      distM = measurable_distance,
      minClusterSize = 1,
      method = "hybrid",
      deepSplit = deep_split,
      verbose = 0
    )))
    # 0 = left unassigned by the tree cut: each such peptide is on its own
    unassigned <- tree_cluster == 0
    tree_cluster[unassigned] <- max(tree_cluster) + seq_len(sum(unassigned))
    cluster[measurable] <- tree_cluster
  } else {
    cluster[measurable] <- stats::cutree(tree, h = cut_height)
  }
  cluster
}

#' Between-peptide variance of profiles, estimated by median matching
#'
#' Peptides of one form do not follow a common profile exactly: there is real
#' peptide-to-peptide spread on top of replicate noise. Its variance is
#' estimated from the residuals of all peptides about their reference
#' profiles, as the value at which a quantile (by default the median) of the
#' absolute standardized residuals equals its expectation for a standard
#' normal. Using the median makes the estimate tolerate a large share of
#' discordant peptides.
#'
#' @param residuals Numeric matrix, peptides x conditions: profile minus
#'   reference.
#' @param noise_variance Numeric matrix, the modelled variance of
#'   \code{residuals} without the between-peptide variance (replicate noise
#'   plus the uncertainty of the reference).
#' @param matching_quantile Quantile to match (0.5 = median).
#' @param verbose Logical.
#'
#' @return List with \code{between_peptide_variance} (0 if replicate noise
#'   already explains the residuals, or with fewer than 20 residuals) and
#'   \code{residuals_used}.
estimate_dispersion <- function(
  residuals,
  noise_variance,
  matching_quantile = 0.5,
  verbose = FALSE
) {
  residual <- as.vector(residuals)
  variance <- as.vector(noise_variance)
  usable <- is.finite(residual) & is.finite(variance) & variance > 0
  residual <- residual[usable]
  variance <- variance[usable]
  if (length(residual) < 20) {
    return(list(between_peptide_variance = 0, residuals_used = length(residual)))
  }

  expected_quantile <- stats::qnorm(1 - (1 - matching_quantile) / 2)
  # Observed minus expected quantile, given a candidate between-peptide variance
  quantile_excess <- function(candidate_variance) {
    stats::quantile(
      abs(residual) / sqrt(variance + candidate_variance),
      matching_quantile,
      names = FALSE
    ) - expected_quantile
  }

  between_peptide_variance <- if (quantile_excess(0) <= 0) {
    # Replicate noise already explains the observed spread
    0
  } else {
    upper_bound <- max(stats::var(residual), 1e-8) * 50
    if (quantile_excess(upper_bound) > 0) {
      upper_bound
    } else {
      tryCatch(
        stats::uniroot(quantile_excess, c(0, upper_bound), tol = 1e-12)$root,
        error = function(e) 0
      )
    }
  }
  if (verbose) {
    message(sprintf(
      "between-peptide variance = %.6f (sd %.4f log2) from %d residuals, matching quantile = %.2f",
      between_peptide_variance,
      sqrt(between_peptide_variance),
      length(residual),
      matching_quantile
    ))
  }
  list(
    between_peptide_variance = between_peptide_variance,
    residuals_used = length(residual)
  )
}

#' Largest mutually consistent set of peptides in each assembly (coherent core)
#'
#' Within each assembly the peptides are partitioned by profile (average
#' linkage on the standardized profile distance, cut at
#' \code{core_cut_height}); the largest cluster is the core. An assembly whose
#' peptides have no measurable distances is all core.
#'
#' @param centered_profiles Numeric matrix, peptides x conditions.
#' @param profile_variance Numeric matrix, peptides x conditions: the variance
#'   the distance is standardized by.
#' @param assembly Character vector, the assembly of each peptide.
#' @param core_cut_height Numeric. Cut height of the partition.
#' @param min_shared_conditions Integer. Conditions two peptides must share for
#'   their distance to be measured.
#'
#' @return List with \code{in_core} (logical, one per peptide) and
#'   \code{n_clusters} (named integer, clusters found per assembly).
find_coherent_core <- function(
  centered_profiles,
  profile_variance,
  assembly,
  core_cut_height = 1.5,
  min_shared_conditions = 2
) {
  in_core <- rep(FALSE, nrow(centered_profiles))
  n_clusters <- integer(0)

  for (assembly_name in unique(assembly)) {
    members <- which(assembly == assembly_name)
    cluster <- if (length(members) < 2) {
      rep(NA_integer_, length(members))
    } else {
      .partition_profiles(
        centered_profiles[members, , drop = FALSE],
        profile_variance[members, , drop = FALSE],
        core_cut_height,
        min_shared_conditions
      )
    }
    if (all(is.na(cluster))) {
      # No measurable structure: every peptide is the reference
      in_core[members] <- TRUE
      n_clusters[assembly_name] <- 1L
      next
    }
    cluster_sizes <- table(cluster)
    in_core[members] <- !is.na(cluster) &
      cluster == as.integer(names(cluster_sizes)[which.max(cluster_sizes)])
    n_clusters[assembly_name] <- length(cluster_sizes)
  }
  list(in_core = in_core, n_clusters = n_clusters)
}

#' Leave-one-out reference profile and offset from the coherent core
#'
#' For each peptide, the mean centered profile and the mean offset of the core
#' members of its assembly other than itself. A peptide that is the core's
#' only member keeps itself in.
#'
#' @param centered_profiles Numeric matrix, peptides x conditions.
#' @param offset Numeric vector, the offset of each peptide.
#' @param in_core Logical vector: the peptide is in its assembly's core.
#' @param assembly Character vector, the assembly of each peptide.
#'
#' @return List with \code{profile} (matrix, the reference profile of each
#'   peptide), \code{offset} (its reference offset), \code{members_averaged}
#'   (matrix: number of core members the reference averages at each condition)
#'   and \code{core_size} (size of the peptide's assembly's core).
core_reference <- function(centered_profiles, offset, in_core, assembly) {
  n_peptides <- nrow(centered_profiles)
  n_conditions <- ncol(centered_profiles)
  reference_profile <- matrix(NA_real_, n_peptides, n_conditions)
  reference_offset <- rep(NA_real_, n_peptides)
  members_averaged <- matrix(0L, n_peptides, n_conditions)
  core_size <- rep(NA_integer_, n_peptides)

  for (assembly_name in unique(assembly)) {
    members <- which(assembly == assembly_name)
    core_members <- members[in_core[members]]
    if (length(core_members) < 1) next
    core_size[members] <- length(core_members)

    core_profiles <- centered_profiles[core_members, , drop = FALSE]
    observed <- is.finite(core_profiles)
    profile_sum <- colSums(ifelse(observed, core_profiles, 0))
    n_observed <- colSums(observed)
    offset_observed <- is.finite(offset[core_members])
    offset_sum <- sum(offset[core_members][offset_observed])
    n_offsets <- sum(offset_observed)

    # Peptides outside the core: mean over the whole core
    reference_profile[members, ] <- matrix(
      profile_sum / n_observed, length(members), n_conditions, byrow = TRUE
    )
    members_averaged[members, ] <- matrix(
      n_observed, length(members), n_conditions, byrow = TRUE
    )
    reference_offset[members] <- offset_sum / n_offsets

    # Core members: leave their own profile out
    if (length(core_members) > 1) {
      own_profile <- ifelse(observed, core_profiles, 0)
      reference_profile[core_members, ] <-
        (matrix(profile_sum, length(core_members), n_conditions, byrow = TRUE) -
          own_profile) /
        (matrix(n_observed, length(core_members), n_conditions, byrow = TRUE) -
          observed)
      members_averaged[core_members, ] <-
        matrix(n_observed, length(core_members), n_conditions, byrow = TRUE) -
        observed
      reference_offset[core_members] <-
        (offset_sum - ifelse(offset_observed, offset[core_members], 0)) /
        (n_offsets - offset_observed)
    }
  }
  reference_profile[!is.finite(reference_profile)] <- NA_real_
  reference_offset[!is.finite(reference_offset)] <- NA_real_
  list(
    profile = reference_profile,
    offset = reference_offset,
    members_averaged = members_averaged,
    core_size = core_size
  )
}

#' Multiple-testing correction of the profile p-values
#'
#' \describe{
#'   \item{\code{"two_stage"}}{Bonferroni within each assembly (capped at 1),
#'     then Benjamini-Hochberg across all assemblies.}
#'   \item{\code{"global_bh"}}{Benjamini-Hochberg across all tested peptides
#'     of all assemblies, with no within-assembly step.}
#' }
#'
#' @param p_values Numeric vector of p-values (\code{NA} for untested
#'   peptides).
#' @param assembly Character vector, the assembly of each peptide.
#' @param alpha Numeric. Threshold on the adjusted p-value.
#' @param method Character. \code{"two_stage"} or \code{"global_bh"}.
#'
#' @return List with \code{within_assembly_p} (the p-values after the
#'   within-assembly step; the raw p-values for \code{"global_bh"}),
#'   \code{adjusted_p} and \code{discordant} (logical: adjusted p-value below
#'   \code{alpha}).
correct_pvalues <- function(
  p_values,
  assembly,
  alpha = 0.05,
  method = c("two_stage", "global_bh")
) {
  method <- match.arg(method)
  tested <- is.finite(p_values)
  tests_in_assembly <- as.numeric(table(assembly[tested])[assembly])
  tests_in_assembly[is.na(tests_in_assembly)] <- 0

  within_assembly_p <- rep(NA_real_, length(p_values))
  within_assembly_p[tested] <- if (method == "two_stage") {
    pmin(p_values[tested] * tests_in_assembly[tested], 1)
  } else {
    p_values[tested]
  }

  adjusted_p <- rep(NA_real_, length(p_values))
  if (any(tested)) {
    adjusted_p[tested] <- stats::p.adjust(within_assembly_p[tested], method = "BH")
  }

  list(
    within_assembly_p = within_assembly_p,
    adjusted_p = adjusted_p,
    discordant = tested & is.finite(adjusted_p) & adjusted_p < alpha
  )
}

#' Group the discordant peptides of each assembly into dCFs
#'
#' Within an assembly, discordant peptides that share a profile form one
#' differential complexoform (dCF). The discordant peptides are clustered
#' hierarchically (Ward.D2) on the standardized profile distance and the tree
#' is cut with the dynamic tree cut (\code{dynamicTreeCut::cutreeDynamic()},
#' hybrid method). The distance is standardized by replicate noise plus the
#' between-peptide variance, so two peptides are one unit apart when they
#' differ by as much as two peptides of the same profile do.
#'
#' The minimum cluster size of the tree cut is one. Clusters of at least
#' \code{min_cluster_size} peptides are labelled \code{dCF1}, \code{dCF2}, ...
#' by decreasing size; smaller ones (a discordant peptide on its own) get
#' \code{singleton_label}; every peptide that is not discordant gets
#' \code{canonical_label}.
#'
#' @param centered_profiles Numeric matrix, peptides x conditions.
#' @param squared_standard_error Numeric matrix, peptides x conditions:
#'   squared standard error of the condition means.
#' @param between_peptide_variance Numeric. Between-peptide variance of
#'   profiles.
#' @param assembly Character vector, the assembly of each peptide.
#' @param discordant Logical vector: the peptide is discordant.
#' @param cut_height Numeric. Height at which two discordant peptides are
#'   still one group when they are the only two (the dynamic tree cut needs
#'   three).
#' @param deep_split Integer 0-4. Sensitivity of the dynamic tree cut to
#'   splitting clusters.
#' @param min_cluster_size Integer. Minimum number of peptides in a labelled
#'   dCF.
#' @param canonical_label Character. Label of peptides that are not
#'   discordant.
#' @param singleton_label Character. Label of discordant peptides outside any
#'   dCF.
#'
#' @return Character vector of labels, one per peptide.
assign_dcf_labels <- function(
  centered_profiles,
  squared_standard_error,
  between_peptide_variance,
  assembly,
  discordant,
  cut_height = 1.5,
  deep_split = 2,
  min_cluster_size = 2,
  canonical_label = "dCF0",
  singleton_label = "dCF-1"
) {
  labels <- rep(canonical_label, nrow(centered_profiles))
  labels[discordant] <- singleton_label

  for (assembly_name in unique(assembly[discordant])) {
    discordant_members <- which(assembly == assembly_name & discordant)
    if (length(discordant_members) < min_cluster_size) next

    cluster <- .partition_profiles(
      centered_profiles[discordant_members, , drop = FALSE],
      squared_standard_error[discordant_members, , drop = FALSE] +
        between_peptide_variance,
      cut_height,
      dynamic = TRUE,
      deep_split = deep_split
    )
    cluster_sizes <- sort(table(cluster), decreasing = TRUE)
    cluster_sizes <- cluster_sizes[cluster_sizes >= min_cluster_size]
    for (size_rank in seq_along(cluster_sizes)) {
      in_cluster <- !is.na(cluster) &
        cluster == as.integer(names(cluster_sizes)[size_rank])
      labels[discordant_members[in_cluster]] <- paste0("dCF", size_rank)
    }
  }
  labels
}

#' Discover complexoforms: discordant peptides and their dCF groups
#'
#' Tests every peptide's profile across conditions against the profile of its
#' assembly's coherent core, and groups the discordant peptides of each
#' assembly into differential complexoforms (dCFs).
#'
#' Steps, for the whole table at once:
#' \enumerate{
#'   \item Completeness filter (\code{prepare_peptides()}).
#'   \item Per-condition means; replicate variance pooled across conditions
#'     and moderated, with one prior per \code{datatype_column} level.
#'   \item Each peptide is centered on its own mean: the offset is set aside,
#'     the profile is tested.
#'   \item Coherent core per assembly; leave-one-out reference profile.
#'   \item Between-peptide variance from all residuals (median matching).
#'   \item Profile test: the squared residuals, each divided by its variance,
#'     are summed and compared with a chi-square on (conditions - 1) degrees
#'     of freedom.
#'   \item Multiple-testing correction.
#'   \item dCF labels for the discordant peptides (dynamic tree cut).
#' }
#'
#' @param data data.table (or coercible) with one row per assembly and peptide,
#'   across all assemblies.
#' @param intensity_columns Character vector of intensity column names (log2).
#' @param assembly_column Character. Column identifying the assembly (the
#'   complex).
#' @param peptide_column Character. Column with peptide identifiers, unique
#'   within an assembly.
#' @param condition_regex Regex with one capture group extracting the condition
#'   label from an intensity column name.
#' @param alpha Numeric. Threshold on the adjusted p-value.
#' @param min_fraction_observed Numeric in [0, 1]. Minimum fraction of
#'   non-missing intensities required to test a peptide.
#' @param min_conditions Integer. Minimum number of conditions that must meet
#'   \code{min_replicates_per_condition}.
#' @param min_replicates_per_condition Integer. Minimum non-missing replicates
#'   for a condition to count.
#' @param datatype_column Character or \code{NULL}. Column to stratify the
#'   variance prior by (e.g. PTM vs non-modified peptides).
#' @param core_cut_height Numeric. Cut height for the coherent core, in units
#'   of the standardized profile distance.
#' @param matching_quantile Numeric. Quantile matched when estimating the
#'   between-peptide variance.
#' @param variance_aware_core Logical. If \code{TRUE}, the core is found a
#'   second time with the between-peptide variance added to the variance the
#'   distance is standardized by, and the reference, the between-peptide
#'   variance and the test are recomputed from that core. With \code{FALSE}
#'   the core distance is in units of replicate noise only.
#' @param correction Character. \code{"two_stage"} (Bonferroni within
#'   assembly, then Benjamini-Hochberg across all assemblies) or
#'   \code{"global_bh"} (Benjamini-Hochberg across all peptides only); see
#'   \code{correct_pvalues()}.
#' @param deep_split Integer 0-4. Sensitivity of the dynamic tree cut that
#'   groups the discordant peptides into dCFs; higher gives more, smaller
#'   groups. See \code{assign_dcf_labels()}.
#' @param min_cluster_size Integer. Minimum number of discordant peptides that
#'   form a labelled dCF. A smaller group - a discordant peptide on its own -
#'   gets \code{singleton_label}.
#' @param canonical_label Character. Label for peptides that are not
#'   discordant.
#' @param singleton_label Character. Label for discordant peptides outside any
#'   dCF group.
#' @param verbose Logical.
#'
#' @return A list with:
#'
#'   - `peptide_results`: data.table with every column and row of `data`, plus
#'     `n_cond_tested` (conditions the test used), `core_size`, `in_core`,
#'     `p_profile` (p-value of the profile test), `p_offset` (p-value of the
#'     peptide's offset against the core's; reported, not used), `p_stage1`
#'     (after the within-assembly step), `p_adj`, `discordant` and `dCF`.
#'     Peptides that fail the filter keep their rows, with `NA` statistics,
#'     `discordant = FALSE` and the canonical label.
#'   - `prior_table`: the fitted variance prior per stratum.
#'   - `tau2`: the between-peptide variance.
#'   - `model`: the fitted noise model, row-aligned with `data`: `tested`
#'     (logical), `cond_mean` (condition means) and `se2` (their squared
#'     standard errors), both peptides x conditions and `NA` for untested
#'     peptides, `s2` (moderated replicate variance per peptide) and
#'     `cond_levels` (the conditions in order).
#'
#' @export
discover_complexoforms <- function(
  data,
  intensity_columns,
  assembly_column = "complex_name",
  peptide_column = "Peptide",
  condition_regex = "^(C_[0-9]+)_R_[0-9]+$",
  alpha = 0.05,
  min_fraction_observed = 0,
  min_conditions = 2,
  min_replicates_per_condition = 1,
  datatype_column = NULL,
  core_cut_height = 1.5,
  matching_quantile = 0.5,
  variance_aware_core = FALSE,
  correction = c("two_stage", "global_bh"),
  deep_split = 2,
  min_cluster_size = 2,
  canonical_label = "dCF0",
  singleton_label = "dCF-1",
  verbose = FALSE
) {
  correction <- match.arg(correction)

  prepared <- prepare_peptides(
    data,
    intensity_columns,
    assembly_column,
    peptide_column,
    condition_regex,
    min_fraction_observed,
    min_conditions,
    min_replicates_per_condition,
    datatype_column,
    verbose
  )

  peptide_results <- data.table::copy(as.data.table(data))
  peptide_results[, n_cond_tested := NA_integer_]
  peptide_results[, core_size := NA_integer_]
  peptide_results[, in_core := NA]
  for (p_value_column in c("p_profile", "p_offset", "p_stage1", "p_adj")) {
    peptide_results[, (p_value_column) := NA_real_]
  }
  peptide_results[, discordant := FALSE]
  peptide_results[, dCF := canonical_label]

  n_rows <- nrow(peptide_results)
  n_conditions <- length(prepared$condition_levels)
  empty_condition_matrix <- matrix(
    NA_real_, n_rows, n_conditions,
    dimnames = list(NULL, prepared$condition_levels)
  )
  model <- list(
    tested = prepared$tested,
    cond_mean = empty_condition_matrix,
    se2 = empty_condition_matrix,
    s2 = rep(NA_real_, n_rows),
    cond_levels = prepared$condition_levels
  )

  if (sum(prepared$tested) < 2) {
    if (verbose) {
      message("Insufficient peptides passed filtering. Returning all as canonical.")
    }
    return(list(
      peptide_results = peptide_results,
      prior_table = data.table(),
      tau2 = NA_real_,
      model = model
    ))
  }

  # ---- noise model ----
  stats_by_condition <- condition_stats(
    prepared$intensities,
    prepared$column_conditions,
    prepared$condition_levels,
    min_replicates_per_condition
  )
  variance <- moderate_variance(
    stats_by_condition$variance,
    stats_by_condition$replicates,
    prior_stratum = prepared$datatype,
    verbose = verbose
  )
  centered <- split_offset_profile(stats_by_condition$mean)

  # Squared standard error of each condition mean, from the moderated variance
  squared_standard_error <- sweep(
    1 / pmax(stats_by_condition$replicates, 1),
    1,
    variance$moderated_variance,
    "*"
  )
  squared_standard_error[!is.finite(stats_by_condition$mean)] <- NA_real_

  # ---- coherent core, reference, between-peptide variance ----
  # The reference is built without a between-peptide variance; its residuals
  # then give that variance.
  fit_reference <- function(in_core) {
    reference <- core_reference(
      centered$profile,
      centered$offset,
      in_core,
      prepared$assembly
    )
    # Replicate noise plus the uncertainty of the reference
    noise_variance <- squared_standard_error *
      (1 + 1 / pmax(reference$members_averaged, 1))
    dispersion <- estimate_dispersion(
      centered$profile - reference$profile,
      noise_variance,
      matching_quantile = matching_quantile,
      verbose = verbose
    )
    list(
      in_core = in_core,
      reference = reference,
      noise_variance = noise_variance,
      between_peptide_variance = dispersion$between_peptide_variance
    )
  }

  reference_fit <- fit_reference(
    find_coherent_core(
      centered$profile,
      squared_standard_error,
      prepared$assembly,
      core_cut_height = core_cut_height
    )$in_core
  )
  if (variance_aware_core) {
    reference_fit <- fit_reference(
      find_coherent_core(
        centered$profile,
        squared_standard_error + reference_fit$between_peptide_variance,
        prepared$assembly,
        core_cut_height = core_cut_height
      )$in_core
    )
  }
  reference <- reference_fit$reference
  between_peptide_variance <- reference_fit$between_peptide_variance

  # ---- profile test: (conditions - 1) degrees of freedom ----
  squared_standardized_residual <- (centered$profile - reference$profile)^2 /
    (reference_fit$noise_variance + between_peptide_variance)
  squared_standardized_residual[!is.finite(squared_standardized_residual)] <- NA_real_
  conditions_tested <- rowSums(is.finite(squared_standardized_residual))
  # One degree of freedom is spent centering each peptide on its own mean
  profile_p <- ifelse(
    conditions_tested >= min_conditions,
    stats::pchisq(
      rowSums(squared_standardized_residual, na.rm = TRUE),
      df = pmax(conditions_tested - 1, 1),
      lower.tail = FALSE
    ),
    NA_real_
  )

  # ---- offset contrast: 1 degree of freedom, reported only ----
  values_observed <- rowSums(
    stats_by_condition$replicates * is.finite(stats_by_condition$mean)
  )
  offset_variance <- variance$moderated_variance / pmax(values_observed, 1) +
    between_peptide_variance / pmax(centered$conditions_observed, 1)
  offset_z <- (centered$offset - reference$offset) / sqrt(offset_variance)
  offset_p <- ifelse(is.finite(offset_z), 2 * stats::pnorm(-abs(offset_z)), NA_real_)

  corrected <- correct_pvalues(
    profile_p,
    prepared$assembly,
    alpha = alpha,
    method = correction
  )

  if (verbose) {
    message(sum(corrected$discordant), " peptides marked as discordant")
  }

  # ---- dCF labels ----
  dcf_labels <- assign_dcf_labels(
    centered$profile,
    squared_standard_error,
    between_peptide_variance,
    prepared$assembly,
    corrected$discordant,
    cut_height = core_cut_height,
    deep_split = deep_split,
    min_cluster_size = min_cluster_size,
    canonical_label = canonical_label,
    singleton_label = singleton_label
  )

  tested_rows <- which(prepared$tested)
  set(peptide_results, tested_rows, "n_cond_tested", as.integer(conditions_tested))
  set(peptide_results, tested_rows, "core_size", as.integer(reference$core_size))
  set(peptide_results, tested_rows, "in_core", reference_fit$in_core)
  set(peptide_results, tested_rows, "p_profile", profile_p)
  set(peptide_results, tested_rows, "p_offset", offset_p)
  set(peptide_results, tested_rows, "p_stage1", corrected$within_assembly_p)
  set(peptide_results, tested_rows, "p_adj", corrected$adjusted_p)
  set(peptide_results, tested_rows, "discordant", corrected$discordant)
  set(peptide_results, tested_rows, "dCF", dcf_labels)

  model$cond_mean[tested_rows, ] <- stats_by_condition$mean
  model$se2[tested_rows, ] <- squared_standard_error
  model$s2[tested_rows] <- variance$moderated_variance

  list(
    peptide_results = peptide_results[],
    prior_table = variance$prior_table,
    tau2 = between_peptide_variance,
    model = model
  )
}

#' Find Number of Clusters from dCF Labels
#'
#' Counts the unique dCF labels, excluding the singleton label: the canonical
#' group plus one per dCF.
#'
#' @param data data.table containing a "dCF" column with dCF labels.
#' @param singleton_label Character. Label that is not counted.
#'
#' @return Integer number of clusters.
#'
#' @export
find_nclust <- function(data, singleton_label = "dCF-1") {
  if (!"dCF" %in% colnames(data)) {
    stop("Input data must contain a 'dCF' column.")
  }
  length(unique(data$dCF[data$dCF != singleton_label]))
}
