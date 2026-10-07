library(data.table)

# Uses parse_conditions() and condition_means() from 04-Clustering/VsClust.R

#' Parse peptide start and end positions
#'
#' Reads the residue range from position strings of the form
#' \code{"ACC [start-end]"}, or a bare \code{"start-end"}. A string that lists
#' several proteins (\code{"ACC1 [s-e]; ACC2 [s-e]"}) is read at the entry of
#' the row's own accession; if that entry is missing the position is
#' \code{NA}.
#'
#' @param positions Character vector of position strings.
#' @param accessions Character vector of accessions, one per position string.
#'
#' @return data.table with integer columns \code{Start} and \code{End}.
parse_peptide_positions <- function(positions, accessions) {
  positions <- as.character(positions)
  own_entry <- positions

  lists_several_proteins <- which(grepl(";", positions, fixed = TRUE))
  for (row in lists_several_proteins) {
    entries <- trimws(strsplit(positions[row], ";", fixed = TRUE)[[1]])
    matching_entries <- entries[startsWith(entries, paste0(accessions[row], " ["))]
    own_entry[row] <- if (length(matching_entries) >= 1L) {
      matching_entries[1L]
    } else {
      NA_character_
    }
  }

  # Regex match of each entry: the whole match, the start and the end
  range_match <- regmatches(
    own_entry,
    regexec("\\[(\\d+)\\s*[-–]\\s*(\\d+)\\]", own_entry, perl = TRUE)
  )
  # Positions given as a bare range ("289-298")
  is_bare_range <- lengths(range_match) != 3 & !is.na(own_entry)
  range_match[is_bare_range] <- regmatches(
    own_entry[is_bare_range],
    regexec("^\\s*(\\d+)\\s*[-–]\\s*(\\d+)\\s*$", own_entry[is_bare_range], perl = TRUE)
  )
  matched_integer <- function(match_index) {
    vapply(
      range_match,
      function(match) {
        if (length(match) == 3) as.integer(match[match_index]) else NA_integer_
      },
      integer(1)
    )
  }
  data.table(Start = matched_integer(2), End = matched_integer(3))
}

#' Create peptide cannot-link constraints
#'
#' Constructs a peptide x peptide cannot-link matrix: two peptides may not
#' share a cluster when they cover the same stretch of the same protein but
#' their profiles across conditions differ. Overlapping peptides that behave
#' differently are different forms of that stretch (e.g. modified and
#' unmodified); overlapping peptides that behave alike (say the singly, doubly
#' and triply phosphorylated versions of a stretch) stay free to share a
#' cluster.
#'
#' A pair is a candidate when both peptides come from the same accession, they
#' overlap by more than \code{min_overlap} of the shorter one, and they share
#' at least \code{min_shared_conditions} conditions.
#'
#' \strong{The test.} The two profiles are centered on the conditions they
#' share. In every shared condition the squared difference between them is
#' divided by its variance (the two squared standard errors plus twice the
#' between-peptide variance), and these are summed. The sum is compared with a
#' chi-square on (shared conditions - 1) degrees of freedom, noncentral by
#' what a difference of exactly \code{tolerance} log2 in every condition would
#' add. So only a difference larger than the tolerance counts, and the size of
#' a change counts, not only its shape.
#'
#' \strong{Two ways to apply it}, chosen with \code{grouping}:
#'
#' \describe{
#'   \item{\code{grouping = TRUE}: pattern groups per position.}{The candidate
#'     pairs tie overlapping peptides into stretches; within a stretch the
#'     peptides are grouped by profile, and a candidate pair is linked when
#'     its two peptides end up in different groups. The test compares the
#'     profiles of whole groups, which are more precise than one peptide's,
#'     and a pair within a group is never linked. See
#'     \code{.link_pattern_groups()}.}
#'   \item{\code{grouping = FALSE}: pair by pair.}{Every candidate pair is
#'     tested on its own and linked when its Benjamini-Hochberg-adjusted
#'     p-value (across the candidate pairs) is below \code{link_threshold}.
#'     The \code{pairs} attribute then has a column \code{zone}:
#'     \code{different} (linked), \code{similar} (the same test in the other
#'     direction: the profiles demonstrably differ by less than the tolerance)
#'     or \code{unclear}.}
#' }
#'
#' Either way a pair is linked only on evidence of a difference: where the
#' data cannot tell, the peptides stay free to share a cluster.
#'
#' @param data A `data.table` (or coercible) with one row per peptide of one
#'   complex.
#' @param id_columns Character vector. Columns pasted together (with `"_"`) to
#'   form the peptide identifier; must match those passed to
#'   `vsclust_on_complex()`.
#' @param intensity_columns Character vector of intensity columns.
#' @param condition_regex Regex with one capture group for the condition label
#'   of an intensity column.
#' @param grouping Logical. `TRUE`: pattern groups per position. `FALSE`: each
#'   candidate pair on its own.
#' @param tolerance Numeric. The difference between two profiles (log2, root
#'   mean square over the shared conditions) that does not count as a
#'   difference. Defaults to the square root of the between-peptide variance
#'   in `noise_model`; 0 links on any difference beyond noise.
#' @param link_threshold Numeric. Adjusted p-value below which peptides are
#'   separated. Default 0.05. With `method` it is the threshold of that rule
#'   (default 0.5 for pearson and ccc, 1 for zrmsd).
#' @param min_shared_conditions Integer. Minimum number of conditions with data
#'   in both peptides. Defaults to a third of the conditions.
#' @param min_overlap Numeric. Minimum overlap, as a fraction of the shorter
#'   peptide, that must be exceeded.
#' @param accession_column Character. Column with the protein accession.
#' @param position_column Character. Column with the position string.
#' @param noise_model List with `squared_standard_error` (numeric matrix, rows
#'   aligned with `data`, one column per condition named by its label: squared
#'   standard error of the condition means) and `between_peptide_variance`
#'   (numeric). These are `model$se2` and `tau2` of
#'   `discover_complexoforms()`.
#' @param method Optional. An alternative rule that replaces the test:
#'   `"pearson"` (Pearson correlation of the condition means below
#'   `link_threshold`), `"zrmsd"` (root mean square difference between
#'   robust-z transformed condition means above it) or `"ccc"` (Lin's
#'   concordance correlation coefficient below it; not invariant to peptide
#'   offsets). These need no `noise_model`. Also accepted: `"groups"`
#'   (`grouping = TRUE`), `"tolerance"` (`grouping = FALSE`) and
#'   `"calibrated"` (`grouping = FALSE`, `tolerance = 0`).
#' @param verbose Logical. Whether to emit progress messages.
#'
#' @return A logical matrix (`TRUE` = cannot-link) with dimnames set to the
#'   peptide identifiers. The attribute `pairs` holds the candidate pairs
#'   (`peptide_i`, `peptide_j`, the number of `shared` conditions, `linked`)
#'   with their `score`: pair by pair the p-value of the test, with grouping
#'   the adjusted p-value at which the two peptides are separated. With
#'   grouping the attribute `groups` gives each peptide's stretch and pattern
#'   group.
#'
#' @export
create_cannot_link <- function(
  data,
  id_columns = c("Gene name", "Peptide"),
  intensity_columns,
  condition_regex,
  grouping = TRUE,
  tolerance = NULL,
  link_threshold = NULL,
  min_shared_conditions = NULL,
  min_overlap = 0.5,
  accession_column = "Accession",
  position_column = "Position in master protein",
  noise_model = NULL,
  method = NULL,
  verbose = FALSE
) {
  if (is.null(method)) {
    method <- if (isTRUE(grouping)) "groups" else "tolerance"
  } else {
    method <- match.arg(
      method,
      c("pearson", "calibrated", "tolerance", "groups", "zrmsd", "ccc")
    )
  }
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  missing_columns <- setdiff(
    c(id_columns, intensity_columns, accession_column, position_column),
    names(data)
  )
  if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  if (is.null(link_threshold)) {
    link_threshold <- c(
      pearson = 0.5, calibrated = 0.05, tolerance = 0.05, groups = 0.05,
      zrmsd = 1, ccc = 0.5
    )[[method]]
  }
  uses_noise_model <- method %in% c("calibrated", "tolerance", "groups")
  squared_standard_error <- noise_model$squared_standard_error
  between_peptide_variance <- noise_model$between_peptide_variance
  if (uses_noise_model) {
    if (is.null(squared_standard_error) || is.null(between_peptide_variance)) {
      stop(
        "create_cannot_link() needs noise_model with squared_standard_error ",
        "and between_peptide_variance (from discover_complexoforms())."
      )
    }
    if (nrow(squared_standard_error) != nrow(data)) {
      stop("noise_model$squared_standard_error must have one row per row of data.")
    }
    if (is.null(tolerance)) {
      tolerance <- sqrt(between_peptide_variance)
    }
  }

  peptide_ids <- do.call(paste, c(data[, id_columns, with = FALSE], sep = "_"))
  if (anyDuplicated(peptide_ids)) {
    stop("Peptide identifiers built from id_columns are not unique.")
  }

  column_conditions <- parse_conditions(intensity_columns, condition_regex)
  intensities <- as.matrix(data[, intensity_columns, with = FALSE])
  storage.mode(intensities) <- "double"
  condition_mean <- condition_means(intensities, column_conditions)
  if (is.null(min_shared_conditions)) {
    min_shared_conditions <- floor(ncol(condition_mean) / 3)
  }
  if (uses_noise_model) {
    if (!all(colnames(condition_mean) %in% colnames(squared_standard_error))) {
      stop(
        "noise_model$squared_standard_error needs one column per condition, ",
        "named by its label."
      )
    }
    # Same condition order as condition_mean
    squared_standard_error <- squared_standard_error[,
      colnames(condition_mean),
      drop = FALSE
    ]
  }

  accession <- as.character(data[[accession_column]])
  position <- parse_peptide_positions(data[[position_column]], accession)
  has_position <- !is.na(position$Start) & !is.na(position$End) & !is.na(accession)

  if (verbose) {
    message("Parsed ", sum(has_position), " of ", nrow(data), " peptide positions")
  }

  # ---- candidate pairs: same accession, overlapping, enough shared conditions
  # first_row and second_row are the rows of the two peptides in `data`
  candidates <- rbindlist(lapply(unique(accession[has_position]), function(protein) {
    protein_rows <- which(has_position & accession == protein)
    if (length(protein_rows) < 2) {
      return(NULL)
    }
    row_pairs <- utils::combn(protein_rows, 2)
    first_row <- row_pairs[1, ]
    second_row <- row_pairs[2, ]
    overlap_length <- pmin(position$End[first_row], position$End[second_row]) -
      pmax(position$Start[first_row], position$Start[second_row]) + 1
    shorter_length <- pmin(
      position$End[first_row] - position$Start[first_row],
      position$End[second_row] - position$Start[second_row]
    ) + 1
    overlapping <- overlap_length > 0 & overlap_length / shorter_length > min_overlap
    data.table(first_row = first_row[overlapping], second_row = second_row[overlapping])
  }))

  cannot_link <- matrix(
    FALSE,
    nrow = length(peptide_ids),
    ncol = length(peptide_ids),
    dimnames = list(peptide_ids, peptide_ids)
  )
  if (!is.null(candidates) && nrow(candidates) > 0) {
    candidates[,
      shared := rowSums(is.finite(condition_mean[first_row, , drop = FALSE]) &
        is.finite(condition_mean[second_row, , drop = FALSE]))
    ]
    candidates <- candidates[shared >= min_shared_conditions]
  }
  if (is.null(candidates) || nrow(candidates) == 0) {
    attr(cannot_link, "pairs") <- data.table()
    return(cannot_link)
  }

  # ---- score the candidate pairs
  robust_z <- function(values) {
    median_absolute_deviation <- mad(values, constant = 1, na.rm = TRUE)
    if (!is.finite(median_absolute_deviation) || median_absolute_deviation == 0) {
      return(rep(NA_real_, length(values)))
    }
    (values - median(values, na.rm = TRUE)) / median_absolute_deviation
  }

  score_pair <- function(first_row, second_row) {
    first_profile <- condition_mean[first_row, ]
    second_profile <- condition_mean[second_row, ]
    shared <- is.finite(first_profile) & is.finite(second_profile)
    if (sum(shared) < 3) {
      return(NA_real_)
    }
    if (uses_noise_model) {
      # Difference between the profiles, each centered on the shared conditions
      difference <- (first_profile[shared] - mean(first_profile[shared])) -
        (second_profile[shared] - mean(second_profile[shared]))
      difference_variance <- squared_standard_error[first_row, shared] +
        squared_standard_error[second_row, shared] +
        2 * between_peptide_variance
    }
    switch(
      method,
      pearson = stats::cor(first_profile[shared], second_profile[shared]),
      calibrated = stats::pchisq(
        sum(difference^2 / difference_variance),
        df = sum(shared) - 1,
        lower.tail = FALSE
      ),
      tolerance = .p_exceeds_tolerance(
        sum(difference^2 / difference_variance),
        sum(shared) - 1,
        tolerance^2 * sum(1 / difference_variance)
      ),
      zrmsd = {
        first_z <- robust_z(first_profile)
        second_z <- robust_z(second_profile)
        shared <- shared & is.finite(first_z) & is.finite(second_z)
        if (sum(shared) < 3) {
          NA_real_
        } else {
          sqrt(mean((first_z[shared] - second_z[shared])^2))
        }
      },
      ccc = {
        denominator <- stats::var(first_profile[shared]) +
          stats::var(second_profile[shared]) +
          (mean(first_profile[shared]) - mean(second_profile[shared]))^2
        if (!is.finite(denominator) || denominator == 0) {
          NA_real_
        } else {
          2 * stats::cov(first_profile[shared], second_profile[shared]) / denominator
        }
      }
    )
  }

  pattern_groups <- NULL
  if (method == "groups") {
    grouped <- .link_pattern_groups(
      condition_mean,
      squared_standard_error,
      between_peptide_variance,
      candidates,
      tolerance,
      peptide_ids,
      alpha = link_threshold
    )
    candidates[, score := grouped$score]
    pattern_groups <- grouped$groups
  } else {
    candidates[, score := mapply(score_pair, first_row, second_row)]
  }

  candidates[,
    linked := switch(
      method,
      pearson = score < link_threshold,
      calibrated = stats::p.adjust(score, method = "BH") < link_threshold,
      tolerance = stats::p.adjust(score, method = "BH") < link_threshold,
      # Already adjusted, across the branchings of all stretches
      groups = score < link_threshold,
      zrmsd = score > link_threshold,
      ccc = score < link_threshold
    )
  ]
  candidates[is.na(linked), linked := FALSE]
  if (method == "tolerance") {
    # 1 - score is the p-value of the test in the other direction
    candidates[, zone := ifelse(
      linked,
      "different",
      ifelse(
        !is.na(score) & stats::p.adjust(1 - score, method = "BH") < link_threshold,
        "similar",
        "unclear"
      )
    )]
  }

  links <- candidates[linked == TRUE]
  cannot_link[cbind(links$first_row, links$second_row)] <- TRUE
  cannot_link[cbind(links$second_row, links$first_row)] <- TRUE

  if (verbose) {
    message(
      "Created ",
      nrow(links),
      " cannot-link constraints from ",
      nrow(candidates),
      " candidate pairs (",
      method,
      ", threshold ",
      link_threshold,
      ")"
    )
  }

  attr(cannot_link, "pairs") <- data.table(
    peptide_i = peptide_ids[candidates$first_row],
    peptide_j = peptide_ids[candidates$second_row],
    shared = candidates$shared,
    score = candidates$score,
    linked = candidates$linked
  )
  if (method == "tolerance") {
    attr(cannot_link, "pairs")[, zone := candidates$zone]
  }
  if (!is.null(pattern_groups)) {
    attr(cannot_link, "groups") <- pattern_groups
  }
  cannot_link
}

#' P-value for "the profiles differ by more than the tolerance"
#'
#' Upper tail of a chi-square distribution, noncentral when there is a
#' tolerance.
#'
#' @param statistic Numeric. Sum of squared standardized differences.
#' @param degrees_of_freedom Numeric.
#' @param noncentrality Numeric. What a difference of exactly the tolerance
#'   adds to the statistic; 0 for no tolerance.
#'
#' @return The p-value, or \code{NA} if the statistic is not finite or there
#'   are no degrees of freedom.
.p_exceeds_tolerance <- function(statistic, degrees_of_freedom, noncentrality) {
  if (!is.finite(statistic) || degrees_of_freedom < 1) {
    return(NA_real_)
  }
  if (noncentrality > 0) {
    stats::pchisq(
      statistic,
      df = degrees_of_freedom,
      ncp = noncentrality,
      lower.tail = FALSE
    )
  } else {
    stats::pchisq(statistic, df = degrees_of_freedom, lower.tail = FALSE)
  }
}

#' Profile of a group of peptides, with the peptides' own levels removed
#'
#' Fits "condition mean = peptide level + group profile" to the group's
#' condition means by weighted least squares, each value weighted by the
#' inverse of its squared standard error plus the between-peptide variance.
#' Only observed values enter, so the fit stays unbiased when the peptides are
#' observed in different conditions.
#'
#' @param condition_mean Numeric matrix, the group's peptides x conditions.
#' @param squared_standard_error Numeric matrix of the same shape.
#' @param between_peptide_variance Numeric.
#'
#' @return List with \code{profile} (per condition, \code{NA} where no member
#'   is observed) and \code{variance} (the variance of each profile value).
.group_profile <- function(
  condition_mean,
  squared_standard_error,
  between_peptide_variance
) {
  weight <- 1 / (squared_standard_error + between_peptide_variance)
  weight[!is.finite(condition_mean) | !is.finite(weight)] <- 0
  condition_mean[weight == 0] <- 0
  peptide_level <- rowSums(weight * condition_mean) / rowSums(weight)
  peptide_level[!is.finite(peptide_level)] <- 0
  profile <- rep(0, ncol(condition_mean))
  # Alternate between the profile and the peptide levels until they settle
  for (iteration in seq_len(50)) {
    profile <- colSums(weight * (condition_mean - peptide_level)) / colSums(weight)
    profile[!is.finite(profile)] <- 0
    updated_level <- rowSums(weight * sweep(condition_mean, 2, profile)) /
      rowSums(weight)
    updated_level[!is.finite(updated_level)] <- 0
    has_settled <- max(abs(updated_level - peptide_level)) < 1e-9
    peptide_level <- updated_level
    if (has_settled) break
  }
  condition_weight <- colSums(weight)
  list(
    profile = ifelse(condition_weight > 0, profile, NA_real_),
    variance = ifelse(condition_weight > 0, 1 / condition_weight, NA_real_)
  )
}

#' Cannot-links from pattern groups per position
#'
#' Groups the overlapping peptides of each stretch by profile, so that
#' overlapping peptides that behave alike may share a cluster and overlapping
#' peptides that behave differently may not.
#'
#' \enumerate{
#'   \item \strong{Stretches.} The candidate pairs tie the peptides into
#'     connected sets: a stretch of the protein and every peptide version that
#'     covers it.
#'   \item \strong{Tree.} The peptides of a stretch are joined into a tree by
#'     average linkage on their pairwise standardized profile difference (the
#'     statistic of the pairwise test per degree of freedom).
#'   \item \strong{Tests.} At every branching of the tree the two branches are
#'     compared as groups: each branch's profile is the weighted fit of
#'     \code{.group_profile()}, and the difference between the two profiles is
#'     tested against its variance and the \code{tolerance}. A branch of
#'     several peptides has a more precise profile than one peptide.
#'   \item \strong{Groups.} The p-values of all branchings, over all stretches,
#'     are adjusted together (Benjamini-Hochberg). Going down from the top of
#'     each tree, a branching splits the peptides when its adjusted p-value is
#'     below \code{alpha}; below a branching that does not split, nothing is
#'     split. The unsplit branches are the pattern groups.
#' }
#'
#' @param condition_mean Numeric matrix, peptides x conditions.
#' @param squared_standard_error Numeric matrix of the same shape.
#' @param between_peptide_variance Numeric.
#' @param candidates data.table of candidate pairs, with the rows of the two
#'   peptides in \code{first_row} and \code{second_row}.
#' @param tolerance Numeric. Difference between two profiles (log2) that does
#'   not count as a difference.
#' @param peptide_ids Character vector, the identifier of each peptide.
#' @param alpha Numeric. Adjusted p-value below which a branching splits.
#'
#' @return List with \code{score}, one value per candidate pair: the adjusted
#'   p-value at which its two peptides are separated (the largest along the
#'   way from the top of the tree to the branching that separates them), so
#'   the pair is in different groups exactly when the score is below
#'   \code{alpha}; and \code{groups}, a data.table with the \code{stretch} and
#'   \code{pattern_group} of every \code{peptide} in a candidate pair.
.link_pattern_groups <- function(
  condition_mean,
  squared_standard_error,
  between_peptide_variance,
  candidates,
  tolerance,
  peptide_ids,
  alpha = 0.05
) {
  # ---- stretches: connected sets of the candidate pairs
  # Every peptide points to another peptide of its stretch, or to itself when
  # it is the stretch's representative
  points_to <- seq_len(nrow(condition_mean))
  stretch_representative <- function(peptide) {
    while (points_to[peptide] != peptide) peptide <- points_to[peptide]
    peptide
  }
  for (pair in seq_len(nrow(candidates))) {
    first_representative <- stretch_representative(candidates$first_row[pair])
    second_representative <- stretch_representative(candidates$second_row[pair])
    if (first_representative != second_representative) {
      points_to[max(first_representative, second_representative)] <-
        min(first_representative, second_representative)
    }
  }
  candidate_peptides <- sort(unique(c(candidates$first_row, candidates$second_row)))
  representative <- vapply(candidate_peptides, stretch_representative, integer(1))
  stretches <- unname(split(candidate_peptides, representative))

  # Statistic of the pairwise test per degree of freedom
  pair_distance <- function(first_row, second_row) {
    shared <- is.finite(condition_mean[first_row, ]) &
      is.finite(condition_mean[second_row, ])
    if (sum(shared) < 3) return(NA_real_)
    difference <-
      (condition_mean[first_row, shared] - mean(condition_mean[first_row, shared])) -
      (condition_mean[second_row, shared] - mean(condition_mean[second_row, shared]))
    difference_variance <- squared_standard_error[first_row, shared] +
      squared_standard_error[second_row, shared] +
      2 * between_peptide_variance
    sum(difference^2 / difference_variance) / (sum(shared) - 1)
  }
  # P-value for "the two branches differ by more than the tolerance"
  branch_test <- function(left_rows, right_rows) {
    left <- .group_profile(
      condition_mean[left_rows, , drop = FALSE],
      squared_standard_error[left_rows, , drop = FALSE],
      between_peptide_variance
    )
    right <- .group_profile(
      condition_mean[right_rows, , drop = FALSE],
      squared_standard_error[right_rows, , drop = FALSE],
      between_peptide_variance
    )
    shared <- is.finite(left$profile) & is.finite(right$profile)
    if (sum(shared) < 3) return(NA_real_)
    inverse_variance <- 1 / (left$variance[shared] + right$variance[shared])
    difference <- left$profile[shared] - right$profile[shared]
    difference <- difference - sum(inverse_variance * difference) / sum(inverse_variance)
    .p_exceeds_tolerance(
      sum(inverse_variance * difference^2),
      sum(shared) - 1,
      tolerance^2 * sum(inverse_variance)
    )
  }

  # ---- a tree per stretch, and a test per branching
  # `merge` is hclust's: row b joins two branches, a negative entry being a
  # single peptide (by its position in the stretch) and a positive one an
  # earlier branching
  trees <- lapply(stretches, function(stretch_rows) {
    n_peptides <- length(stretch_rows)
    distance <- matrix(0, n_peptides, n_peptides)
    for (first in seq_len(n_peptides - 1)) {
      for (second in (first + 1):n_peptides) {
        distance[first, second] <- distance[second, first] <-
          pair_distance(stretch_rows[first], stretch_rows[second])
      }
    }
    # Pairs that cannot be compared are far apart
    distance[!is.finite(distance)] <- 2 * max(c(distance[is.finite(distance)], 1))
    tree <- stats::hclust(stats::as.dist(distance), method = "average")
    peptides_below <- vector("list", n_peptides - 1)
    branching_p <- rep(NA_real_, n_peptides - 1)
    for (branching in seq_len(n_peptides - 1)) {
      branches <- lapply(tree$merge[branching, ], function(branch) {
        if (branch < 0) -branch else peptides_below[[branch]]
      })
      peptides_below[[branching]] <- c(branches[[1]], branches[[2]])
      branching_p[branching] <- branch_test(
        stretch_rows[branches[[1]]],
        stretch_rows[branches[[2]]]
      )
    }
    list(
      rows = stretch_rows,
      merge = tree$merge,
      peptides_below = peptides_below,
      branching_p = branching_p
    )
  })

  # ---- adjust over all branchings, then split from the top of each tree
  all_branching_p <- unlist(lapply(trees, `[[`, "branching_p"))
  all_branching_p[is.na(all_branching_p)] <- 1
  all_adjusted_p <- stats::p.adjust(all_branching_p, method = "BH")
  branchings_before <- cumsum(c(0, lengths(lapply(trees, `[[`, "branching_p"))))

  pair_score <- rep(NA_real_, nrow(candidates))
  position_in_stretch <- integer(nrow(condition_mean))
  group_tables <- vector("list", length(trees))
  for (stretch in seq_along(trees)) {
    tree <- trees[[stretch]]
    n_peptides <- length(tree$rows)
    top_branching <- n_peptides - 1
    adjusted_p <- all_adjusted_p[branchings_before[stretch] + seq_len(top_branching)]
    # Largest adjusted p-value from the top of the tree down to each branching:
    # a branching splits when this is below alpha
    separation_p <- rep(NA_real_, top_branching)
    separation_p[top_branching] <- adjusted_p[top_branching]
    pattern_group <- integer(n_peptides)
    groups_assigned <- 0L
    for (branching in rev(seq_len(top_branching))) {
      splits <- separation_p[branching] < alpha
      for (branch in tree$merge[branching, ]) {
        if (branch > 0) {
          separation_p[branch] <- max(separation_p[branching], adjusted_p[branch])
        }
        # A single peptide below a split is a group of its own
        if (splits && branch < 0) {
          groups_assigned <- groups_assigned + 1L
          pattern_group[-branch] <- groups_assigned
        }
        # A branch below a split that does not split itself is one group
        if (splits && branch > 0 && !(separation_p[branch] < alpha)) {
          groups_assigned <- groups_assigned + 1L
          pattern_group[tree$peptides_below[[branch]]] <- groups_assigned
        }
      }
    }
    if (!(separation_p[top_branching] < alpha)) pattern_group[] <- 1L

    # The branching that separates two peptides is the first to hold both
    separating_branching <- matrix(NA_integer_, n_peptides, n_peptides)
    for (branching in seq_len(top_branching)) {
      branches <- lapply(tree$merge[branching, ], function(branch) {
        if (branch < 0) -branch else tree$peptides_below[[branch]]
      })
      separating_branching[branches[[1]], branches[[2]]] <- branching
      separating_branching[branches[[2]], branches[[1]]] <- branching
    }
    position_in_stretch[tree$rows] <- seq_len(n_peptides)
    pairs_in_stretch <- which(candidates$first_row %in% tree$rows)
    pair_score[pairs_in_stretch] <- separation_p[separating_branching[cbind(
      position_in_stretch[candidates$first_row[pairs_in_stretch]],
      position_in_stretch[candidates$second_row[pairs_in_stretch]]
    )]]
    group_tables[[stretch]] <- data.table(
      peptide = peptide_ids[tree$rows],
      stretch = stretch,
      pattern_group = pattern_group
    )
  }
  list(score = pair_score, groups = rbindlist(group_tables))
}


#' Create VSClust Restriction Matrix from Initial Clustering Results
#'
#' Runs an unconstrained VSClust clustering, then combines its cluster
#' assignments with the cannot-link constraints into a peptide x cluster
#' restriction matrix for a second, constrained clustering.
#'
#' A restriction names a cluster of the unconstrained clustering. The
#' constrained clustering therefore has to start from the unconstrained
#' clustering's centers, so that its cluster j is the same cluster j: pass
#' both `restriction_matrix` and `initial_centers` from the result to
#' `vsclust_on_complex()`.
#'
#' @param data data.table with intensity values and peptide identifiers
#' @param cannotlink_matrix Peptide x peptide cannot-link constraint matrix
#'   (from `create_cannot_link()`)
#' @param n_clusters Integer. Number of clusters to find
#' @param id_columns Character vector. Columns pasted into peptide identifiers
#' @param intensity_columns Character vector. Columns containing intensity
#'   values
#' @param condition_regex Regex with one capture group for the condition label
#' @param ... Further arguments passed to `vsclust_on_complex()` for the
#'   unconstrained clustering (`standard_deviations`, `scaling`, `n_starts`,
#'   `cores`, `seed`).
#' @param verbose Logical. Whether to print progress messages
#'
#' @return List with three elements:
#'   - restriction_matrix: Matrix with peptides as rows, clusters as columns
#'     * Values are TRUE (forbidden) or FALSE (allowed)
#'     * Rownames = peptide identifiers
#'     * Colnames = cluster numbers (1, 2, 3, ..., n_clusters)
#'   - initial_centers: Centers of the unconstrained clustering (clusters x
#'     conditions), in the cluster order the restrictions refer to
#'   - initial_clustering: The unconstrained VSClust clustering results
#'
#' @details
#' Algorithm:
#' 1. Run unconstrained VSClust clustering with n_clusters
#' 2. Initialize restriction matrix with all FALSE (all peptides allowed everywhere)
#' 3. For each cannot-link pair (pep1, pep2):
#'    a. If assigned to different clusters (cluster1, cluster2):
#'       - Set restriction[pep1, cluster2] = TRUE (pep1 forbidden from cluster2)
#'       - Set restriction[pep2, cluster1] = TRUE (pep2 forbidden from cluster1)
#'    b. If assigned to same cluster (conflict):
#'       - Compare membership values for both peptides in that cluster
#'       - Keep peptide with higher membership, forbid the other
#' 4. A peptide forbidden from every cluster is allowed into all
#'
#' @export
vsclust_to_restrictions <- function(
  data,
  cannotlink_matrix,
  n_clusters,
  id_columns = c("Gene name", "Peptide"),
  intensity_columns,
  condition_regex,
  ...,
  verbose = FALSE
) {
  if (verbose) {
    message("Step 1: Running initial unconstrained VSClust clustering...")
  }

  # Run initial unconstrained clustering
  initial_clustering <- vsclust_on_complex(
    data = data,
    id_columns = id_columns,
    intensity_columns = intensity_columns,
    condition_regex = condition_regex,
    n_clusters = n_clusters,
    restriction_matrix = NULL,
    ...
  )

  # Cluster assignments and memberships, named by peptide identifier. They are
  # NOT in the row order of `data` (runClustWrapper_custom() sorts them), so
  # every lookup below is by name.
  best_clustering <- initial_clustering$ClustOut$Bestcl
  cluster_assignments <- best_clustering$cluster
  membership_matrix <- best_clustering$membership
  peptide_ids <- names(cluster_assignments)

  if (verbose) {
    message(
      "Step 2: Building restriction matrix from clustering + cannot-link..."
    )
    message(
      "  Cluster sizes: ",
      paste(table(cluster_assignments), collapse = ", ")
    )
  }

  # Initialize restriction matrix: all FALSE (all allowed)
  restriction_matrix <- matrix(
    FALSE,
    nrow = length(peptide_ids),
    ncol = n_clusters,
    dimnames = list(peptide_ids, seq_len(n_clusters))
  )

  # Find all cannot-link pairs (upper triangle: each pair once)
  cannotlink_pairs <- which(
    cannotlink_matrix & upper.tri(cannotlink_matrix),
    arr.ind = TRUE
  )

  conflicts_resolved <- 0

  for (pair in seq_len(nrow(cannotlink_pairs))) {
    first_peptide <- rownames(cannotlink_matrix)[cannotlink_pairs[pair, 1]]
    second_peptide <- rownames(cannotlink_matrix)[cannotlink_pairs[pair, 2]]

    # Skip if peptides not in clustering results
    if (!first_peptide %in% peptide_ids || !second_peptide %in% peptide_ids) {
      next
    }

    first_cluster <- cluster_assignments[[first_peptide]]
    second_cluster <- cluster_assignments[[second_peptide]]

    if (first_cluster != second_cluster) {
      # Case A: Peptides in different clusters
      # Forbid each peptide from the other's cluster
      restriction_matrix[first_peptide, second_cluster] <- TRUE
      restriction_matrix[second_peptide, first_cluster] <- TRUE
    } else {
      # Case B: Peptides in same cluster (conflict!)
      # Allow only the peptide with higher membership value
      first_membership <- membership_matrix[first_peptide, first_cluster]
      second_membership <- membership_matrix[second_peptide, first_cluster]
      forbidden_peptide <- if (first_membership > second_membership) {
        second_peptide
      } else {
        first_peptide
      }
      restriction_matrix[forbidden_peptide, first_cluster] <- TRUE
      conflicts_resolved <- conflicts_resolved + 1
    }
  }

  # Handle edge case: peptides forbidden from all clusters
  forbidden_everywhere <- rowSums(restriction_matrix) == n_clusters
  if (any(forbidden_everywhere)) {
    if (verbose) {
      message(
        "  WARNING: ",
        sum(forbidden_everywhere),
        " peptides forbidden from all clusters - allowing into all"
      )
    }
    restriction_matrix[forbidden_everywhere, ] <- FALSE
  }

  # Report statistics
  if (verbose) {
    message(
      "  ",
      nrow(cannotlink_pairs),
      " cannot-link pairs -> ",
      sum(restriction_matrix),
      " restriction entries"
    )
    if (conflicts_resolved > 0) {
      message("  Resolved ", conflicts_resolved, " same-cluster conflicts")
    }
  }

  # Centers of the unconstrained clustering, in the cluster order used above
  initial_centers <- best_clustering$centers
  rownames(initial_centers) <- seq_len(n_clusters)

  list(
    restriction_matrix = restriction_matrix,
    initial_centers = initial_centers,
    initial_clustering = initial_clustering
  )
}

#' Share of cannot-link pairs that share a cluster
#'
#' @param cluster Named integer vector of hard cluster assignments
#'   (`ClustOut$Bestcl$cluster`).
#' @param cannotlink_matrix Peptide x peptide cannot-link matrix.
#'
#' @return Named numeric vector: `pairs` (number of cannot-link pairs with
#'   both peptides clustered) and `violated` (share of them in one cluster).
#'
#' @export
cannot_link_violations <- function(cluster, cannotlink_matrix) {
  linked_pairs <- which(
    cannotlink_matrix & upper.tri(cannotlink_matrix),
    arr.ind = TRUE
  )
  first_cluster <- cluster[rownames(cannotlink_matrix)[linked_pairs[, 1]]]
  second_cluster <- cluster[rownames(cannotlink_matrix)[linked_pairs[, 2]]]
  both_clustered <- !is.na(first_cluster) & !is.na(second_cluster)
  c(
    pairs = sum(both_clustered),
    violated = if (any(both_clustered)) {
      mean(first_cluster[both_clustered] == second_cluster[both_clustered])
    } else {
      NA_real_
    }
  )
}
