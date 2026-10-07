library(data.table)
library(vsclust)
library(limma)
library(matrixStats)
library(parallel)

#' Condition label of each intensity column
#'
#' @param intensity_columns Character vector of intensity column names.
#' @param condition_regex Regex with one capture group for the condition
#'   label.
#' @return Character vector of condition labels, one per column.
parse_conditions <- function(intensity_columns, condition_regex) {
  column_conditions <- sub(condition_regex, "\\1", intensity_columns, perl = TRUE)
  if (any(column_conditions == intensity_columns)) {
    stop(
      "Could not parse conditions from intensity_columns with condition_regex='",
      condition_regex,
      "'"
    )
  }
  column_conditions
}

#' Condition means of an intensity matrix
#'
#' @param intensities Numeric matrix, features x samples.
#' @param column_conditions Character vector, condition label of each column
#'   of \code{intensities}.
#' @return Numeric matrix, features x conditions (in order of first
#'   appearance); \code{NA} where a feature has no value in a condition.
condition_means <- function(intensities, column_conditions) {
  condition_levels <- unique(column_conditions)
  means <- vapply(
    condition_levels,
    function(condition) {
      rowMeans(intensities[, column_conditions == condition, drop = FALSE], na.rm = TRUE)
    },
    numeric(nrow(intensities))
  )
  means <- matrix(
    means,
    nrow = nrow(intensities),
    dimnames = list(rownames(intensities), condition_levels)
  )
  means[is.nan(means)] <- NA_real_
  means
}

#' Prepare Data for VSClust
#'
#' Builds the matrix VSClust clusters: one row per feature with its condition
#' means, followed by the feature's standard deviation in the last column
#' (\code{Sds}), which sets VSClust's variance-sensitive fuzzifier.
#'
#' Conditions are read from the column names with \code{condition_regex}, so
#' the columns can be in any order and the design need not be balanced.
#'
#' @param data data.table with identifier and intensity columns.
#' @param id_columns Character vector of columns pasted with "_" into the row
#'   identifier. The identifiers must be unique.
#' @param intensity_columns Character vector of intensity column names.
#' @param condition_regex Regex with one capture group for the condition
#'   label.
#' @param standard_deviations Optional numeric vector, one standard deviation
#'   per row of \code{data}. If \code{NULL}, the standard deviations are
#'   estimated from \code{data} with limma (moderated, as
#'   \code{vsclust::SignAnalysis()} does).
#'
#' @return Numeric matrix with row identifiers as row names, one column per
#'   condition and a final column \code{Sds}.
#'
#' @export
prepare_for_vsclust <- function(
  data,
  id_columns,
  intensity_columns,
  condition_regex,
  standard_deviations = NULL
) {
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  missing_columns <- setdiff(c(id_columns, intensity_columns), names(data))
  if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }

  row_ids <- do.call(paste, c(data[, id_columns, with = FALSE], sep = "_"))
  if (anyDuplicated(row_ids)) {
    stop("Row identifiers built from id_columns are not unique.")
  }

  intensities <- as.matrix(data[, intensity_columns, with = FALSE])
  storage.mode(intensities) <- "double"
  rownames(intensities) <- row_ids
  column_conditions <- parse_conditions(intensity_columns, condition_regex)

  if (is.null(standard_deviations)) {
    # Moderated standard deviation per feature
    design <- stats::model.matrix(
      ~ 0 + factor(column_conditions, levels = unique(column_conditions))
    )
    variance_fit <- suppressWarnings(
      limma::eBayes(limma::lmFit(intensities, design))
    )
    standard_deviations <- sqrt(variance_fit$s2.post)
  } else if (length(standard_deviations) != nrow(data)) {
    stop("standard_deviations must have one value per row of data.")
  }

  cbind(condition_means(intensities, column_conditions), Sds = standard_deviations)
}

#' Run VSClust on a prepared matrix
#'
#' Scales the condition means, runs VSClust with a fixed number of clusters
#' and returns the clustering in the layout the plotting functions expect.
#'
#' Without \code{initial_centers}, VSClust is started \code{n_starts} times
#' from randomly chosen features and the best solution is kept. With
#' \code{initial_centers}, it is started once, from those centers, so that
#' cluster j of the result grows out of center j.
#'
#' The clusters are numbered by size (largest first) and the features sorted
#' by their highest membership.
#'
#' @param prepared_matrix Numeric matrix from \code{prepare_for_vsclust()}.
#' @param n_clusters Integer. Number of clusters.
#' @param scaling Character. "center" (each feature centered on its own mean,
#'   so the size of its change is kept), "standardize" (centered and scaled
#'   to unit variance: shape only) or "none".
#' @param constraints Optional logical matrix, features x clusters;
#'   \code{TRUE} forbids the feature from the cluster.
#' @param initial_centers Optional numeric matrix, clusters x conditions, in
#'   the scaled space of \code{prepared_matrix}.
#' @param n_starts Integer. Number of random starts.
#' @param cores Integer. Number of worker processes for the random starts.
#' @param seed Optional integer. Makes the random starts reproducible.
#' @param verbose Logical.
#'
#' @return List with \code{dat} (scaled data, rows ordered as \code{Bestcl}),
#'   \code{Bestcl} (VSClust's result: \code{cluster}, \code{membership},
#'   \code{centers}, ...), \code{outFileClust} (the same as \code{dat}),
#'   \code{ClustInd} (members per cluster) and \code{converged} (\code{FALSE}
#'   if VSClust stopped at its limit of 1000 iterations; a warning is given).
#'
#' @export
runClustWrapper_custom <- function(
  prepared_matrix,
  n_clusters,
  scaling = "center",
  constraints = NULL,
  initial_centers = NULL,
  n_starts = 16,
  cores = 1,
  seed = NULL,
  verbose = FALSE
) {
  profiles <- prepared_matrix[, seq_len(ncol(prepared_matrix) - 1), drop = FALSE]
  standard_deviations <- prepared_matrix[, ncol(prepared_matrix)]

  if (!any(scaling == c("standardize", "center", "none"))) {
    stop("parameter scaling needs to be standardize, center or none!")
  }
  # Standardizing divides each feature by its standard deviation across
  # conditions, so its replicate standard deviation is divided by the same
  if ((scaling == "standardize")) {
    standard_deviations <- standard_deviations /
      rowSds(as.matrix(profiles), na.rm = TRUE)
  }
  profiles <- t(scale(
    t(profiles),
    center = (scaling != "none"),
    scale = (scaling == "standardize")
  ))

  # Handle single cluster case - all peptides in one cluster
  if (n_clusters == 1) {
    best_clustering <- list(
      centers = matrix(
        colMeans(profiles, na.rm = TRUE),
        nrow = 1,
        dimnames = list("Cluster 1", colnames(profiles))
      ),
      size = nrow(profiles),
      cluster = setNames(rep(1L, nrow(profiles)), rownames(profiles)),
      membership = matrix(
        1,
        nrow = nrow(profiles),
        ncol = 1,
        dimnames = list(rownames(profiles), "membership of cluster 1")
      )
    )
    return(list(
      dat = profiles,
      Bestcl = best_clustering,
      outFileClust = profiles,
      ClustInd = data.frame(Cluster = "1", Members = nrow(profiles)),
      converged = TRUE
    ))
  }

  if (!is.null(initial_centers)) {
    fuzzifier <- determine_fuzz(dim(profiles), n_clusters, standard_deviations)$m
    best_clustering <- vsclust_algorithm(
      profiles,
      centers = initial_centers,
      m = fuzzifier,
      constraints = constraints,
      iterMax = 1000,
      verbose = verbose
    )
  } else {
    workers <- makeCluster(cores)
    on.exit(stopCluster(workers), add = TRUE)
    if (!is.null(seed)) {
      clusterSetRNGStream(workers, seed)
    }
    clusterExport(
      cl = workers,
      varlist = c("ClustComp", "vsclust_algorithm"),
      envir = environment()
    )
    best_clustering <- ClustComp(
      profiles,
      NClust = n_clusters,
      Sds = standard_deviations,
      constraints = constraints,
      NSs = n_starts,
      cl = workers,
      verbose = verbose
    )$Bestcl
  }

  # VSClust (ClustComp() and the run from given centers above) stops after
  # 1000 iterations whether or not the solution has settled
  converged <- best_clustering$iter < 1000
  if (!converged) {
    warning(
      "VSClust stopped at its limit of 1000 iterations without converging (",
      nrow(profiles),
      " features, ",
      n_clusters,
      " clusters, scaling = \"",
      scaling,
      "\")."
    )
  }

  # Number clusters by size (largest first)
  best_clustering <- SwitchOrder(best_clustering, n_clusters)

  # Sort the features by their highest membership
  membership_order <- order(rowMaxs(best_clustering$membership, na.rm = TRUE))
  best_clustering$cluster <- best_clustering$cluster[membership_order]
  best_clustering$membership <- best_clustering$membership[
    membership_order, ,
    drop = FALSE
  ]
  profiles <- profiles[names(best_clustering$cluster), , drop = FALSE]

  colnames(best_clustering$membership) <-
    paste("membership of cluster", colnames(best_clustering$membership))
  rownames(best_clustering$centers) <-
    paste("Cluster", rownames(best_clustering$centers))
  # Members (membership above 0.5) per cluster
  cluster_sizes <- as.data.frame(table(
    best_clustering$cluster[rowMaxs(best_clustering$membership) > 0.5]
  ))
  if (ncol(cluster_sizes) == 2) {
    colnames(cluster_sizes) <- c("Cluster", "Members")
  } else {
    cluster_sizes <- cbind(
      seq_len(max(best_clustering$cluster)),
      rep(0, max(best_clustering$cluster))
    )
  }

  list(
    dat = profiles,
    Bestcl = best_clustering,
    outFileClust = profiles,
    ClustInd = cluster_sizes,
    converged = converged
  )
}

#' Run VSClust on a single complex-level peptide table
#'
#' Prepares the peptide intensity data of one complex, optionally applies
#' peptide-to-cluster restrictions, and runs VSClust with a fixed number of
#' clusters.
#'
#' @param data A `data.table` (or `data.frame`) containing peptide metadata and
#'   intensity columns for one complex.
#' @param id_columns Character vector of column names used to build a unique
#'   peptide identifier for VSClust.
#' @param intensity_columns Character vector of intensity columns (columns to
#'   be used for clustering).
#' @param condition_regex Regex with one capture group for the condition label
#'   of an intensity column.
#' @param n_clusters Integer. Number of clusters.
#' @param restriction_matrix Optional logical matrix encoding peptide-to-cluster
#'   restrictions, with peptide identifiers as row names. `TRUE` forbids the
#'   peptide from the cluster. Peptides without a row are unrestricted.
#' @param initial_centers Optional numeric matrix of starting centers (clusters
#'   x conditions). With a `restriction_matrix`, pass the centers of the
#'   clustering the restrictions were built against (see
#'   `vsclust_to_restrictions()`), so that the restrictions refer to the same
#'   clusters.
#' @param standard_deviations Optional numeric vector, one standard deviation
#'   per row of `data`; see `prepare_for_vsclust()`.
#' @param scaling Character. "center", "standardize" or "none"; see
#'   `runClustWrapper_custom()`.
#' @param n_starts Integer. Number of random starts when `initial_centers` is
#'   not given.
#' @param cores Integer or `NULL`. Worker processes for the random starts;
#'   defaults to all cores but one.
#' @param seed Optional integer. Makes the random starts reproducible.
#' @param verbose Logical.
#'
#' @return A named list with:
#' \describe{
#'   \item{`ClustOut`}{VSClust output object from
#'   `runClustWrapper_custom()`.}
#'   \item{`original_data`}{Original input data (post coercion to
#'   `data.table`).}
#' }
#'
#' @export
vsclust_on_complex <- function(
  data,
  id_columns,
  intensity_columns,
  condition_regex,
  n_clusters,
  restriction_matrix = NULL,
  initial_centers = NULL,
  standard_deviations = NULL,
  scaling = "center",
  n_starts = 16,
  cores = NULL,
  seed = NULL,
  verbose = FALSE
) {
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  if (is.null(n_clusters)) {
    stop("n_clusters must be provided.")
  }
  if (nrow(data) < 3) {
    stop(
      "Too few rows (",
      nrow(data),
      "). Cannot proceed with clustering."
    )
  }
  if (is.null(cores)) {
    cores <- max(parallel::detectCores() - 1, 1)
  }

  prepared_matrix <- prepare_for_vsclust(
    data,
    id_columns,
    intensity_columns,
    condition_regex,
    standard_deviations
  )

  # Align the restriction matrix to the rows of the prepared matrix by
  # identifier. Peptides without restrictions are allowed in all clusters.
  if (!is.null(restriction_matrix)) {
    restricted_peptides <- intersect(
      rownames(prepared_matrix),
      rownames(restriction_matrix)
    )
    if (length(restricted_peptides) == 0) {
      warning(
        "No matching peptides between data and restriction matrix. Clustering without constraints."
      )
      restriction_matrix <- NULL
    } else {
      aligned_restrictions <- matrix(
        FALSE,
        nrow = nrow(prepared_matrix),
        ncol = ncol(restriction_matrix),
        dimnames = list(rownames(prepared_matrix), colnames(restriction_matrix))
      )
      aligned_restrictions[restricted_peptides, ] <- restriction_matrix[
        restricted_peptides,
      ]
      restriction_matrix <- aligned_restrictions
    }
  }

  if (verbose) {
    message(
      "VSClust: ",
      nrow(prepared_matrix),
      " features, ",
      n_clusters,
      " clusters, ",
      if (is.null(restriction_matrix)) 0 else sum(restriction_matrix),
      " restrictions"
    )
  }

  clustering <- runClustWrapper_custom(
    prepared_matrix,
    n_clusters = n_clusters,
    scaling = scaling,
    constraints = restriction_matrix,
    initial_centers = initial_centers,
    n_starts = n_starts,
    cores = cores,
    seed = seed,
    verbose = FALSE
  )

  list(
    ClustOut = clustering,
    original_data = data
  )
}
