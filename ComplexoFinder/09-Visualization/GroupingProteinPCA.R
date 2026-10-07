library(data.table)
library(ggplot2)
library(cowplot)
library(pcaMethods)
library(scico)

#' PPCA Plot of Peptides Colored by Protein and Cluster Membership
#'
#' Runs probabilistic PCA (PPCA) on condition-averaged peptide intensities and
#' returns two side-by-side scatter plots: one colored by protein of origin
#' (Gene name) and one colored by VSClust cluster assignment (determined by
#' highest membership value) using the scico batlow palette.
#'
#' @param vsclust_result Result from vsclust_on_complex, must contain
#'   \code{$original_data} (data.table with Gene name, Peptide, dCF columns
#'   and intensity columns) and \code{$ClustOut$Bestcl$cluster} and
#'   \code{$ClustOut$Bestcl$membership}.
#' @param intensity_cols Character vector of raw intensity column names.
#' @param n_rep Integer. Number of replicates per condition.
#' @param cond_regex Character. Regex with one capture group for the condition
#'   label (e.g. \code{"^H9_B[0-9]+_(D[0-9]+)$"}).
#' @param complex_name Character. Used as the plot title.
#' @param n_pcs Integer. Number of PCs to compute via PPCA (default 2).
#'
#' @return A cowplot grid with two ggplot panels.
#'
#' @export
plot_grouping_pca <- function(
  vsclust_result,
  intensity_cols,
  n_rep,
  cond_regex,
  complex_name = "",
  n_pcs = 2
) {
  dt <- data.table::copy(vsclust_result$original_data)

  # ── 1. Average replicates per condition ──────────────────────────────────────
  long <- data.table::melt(
    dt,
    id.vars = c("Gene name", "Peptide", "dCF"),
    measure.vars = intensity_cols,
    variable.name = "Sample",
    value.name = "Intensity"
  )

  long[, Condition := sub(cond_regex, "\\1", as.character(Sample))]

  avg <- long[,
    .(Intensity = mean(Intensity, na.rm = TRUE)),
    by = .(`Gene name`, Peptide, dCF, Condition)
  ]

  wide <- data.table::dcast(
    avg,
    `Gene name` + Peptide + dCF ~ Condition,
    value.var = "Intensity"
  )

  # ── 2. Build numeric matrix (rows = peptides, cols = conditions) ─────────────
  cond_cols <- setdiff(names(wide), c("Gene name", "Peptide", "dCF"))
  X <- as.matrix(wide[, cond_cols, with = FALSE])
  # Row IDs match the identifier format used by vsclust_on_complex
  row_ids <- make.unique(paste(wide[["Gene name"]], wide[["Peptide"]], sep = "_"))
  rownames(X) <- row_ids

  # Replace NaN with NA (pcaMethods requires NA, not NaN)
  X[is.nan(X)] <- NA

  # Drop rows with all-NA (cannot impute)
  all_na <- rowSums(!is.na(X)) == 0
  if (any(all_na)) {
    X <- X[!all_na, , drop = FALSE]
    wide <- wide[!all_na]
    row_ids <- row_ids[!all_na]
  }

  # ── 3. PPCA ──────────────────────────────────────────────────────────────────
  pca_res <- pcaMethods::pca(
    X,
    method = "ppca",
    nPcs = n_pcs,
    center = TRUE,
    scale = "uv",
    seed = 42
  )

  scores <- as.data.table(pcaMethods::scores(pca_res), keep.rownames = "row_id")
  var_expl <- pca_res@R2

  # ── 4. Attach Gene name ───────────────────────────────────────────────────────
  scores[, `Gene name` := wide[["Gene name"]][match(row_id, row_ids)]]

  # ── 5. Attach cluster assignment and max membership from constrained VSClust ──
  membership_mat <- vsclust_result$ClustOut$Bestcl$membership
  cluster_vec    <- vsclust_result$ClustOut$Bestcl$cluster

  # max membership value (highest confidence for the assigned cluster)
  max_membership <- matrixStats::rowMaxs(membership_mat, na.rm = TRUE)
  names(max_membership) <- rownames(membership_mat)

  scores[, cluster := cluster_vec[row_id]]
  scores[, max_mem  := max_membership[row_id]]

  # Cluster as factor for discrete coloring; keep numeric order
  n_clusters <- length(unique(stats::na.omit(cluster_vec)))
  cluster_levels <- as.character(seq_len(n_clusters))
  scores[, cluster_fct := factor(as.character(cluster), levels = cluster_levels)]

  # ── 6. Colour palettes ───────────────────────────────────────────────────────
  genes <- sort(unique(scores[["Gene name"]]))
  gene_pal <- scales::hue_pal()(length(genes))
  names(gene_pal) <- genes

  # Batlow discrete palette for clusters
  cluster_pal <- scico::scico(n_clusters, palette = "batlow")
  names(cluster_pal) <- cluster_levels

  pc1_lab <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
  pc2_lab <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])

  base_theme <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  # ── 7. Plot 1: color by Gene name ────────────────────────────────────────────
  p_gene <- ggplot(
    scores,
    aes(x = PC1, y = PC2, color = `Gene name`)
  ) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_color_manual(values = gene_pal) +
    labs(
      x = pc1_lab,
      y = pc2_lab,
      color = "Protein (Gene name)",
      title = "Protein of origin"
    ) +
    coord_fixed() +
    base_theme

  # ── 8. Plot 2: color by cluster (highest membership), batlow palette ──────────
  p_clust <- ggplot(
    scores,
    aes(x = PC1, y = PC2, color = cluster_fct, alpha = max_mem)
  ) +
    geom_point(size = 2.2) +
    scale_color_manual(
      values = cluster_pal,
      na.value = "#CCCCCC",
      drop = FALSE
    ) +
    scale_alpha_continuous(range = c(0.25, 1), limits = c(0, 1)) +
    labs(
      x = pc1_lab,
      y = pc2_lab,
      color = "Cluster",
      alpha = "Max membership",
      title = "Cluster assignment (highest membership)"
    ) +
    coord_fixed() +
    base_theme

  # ── 9. Combine ────────────────────────────────────────────────────────────────
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(
      complex_name,
      fontface = "bold",
      size = 12,
      x = 0.5,
      hjust = 0.5
    )

  combined <- cowplot::plot_grid(p_gene, p_clust, ncol = 2, align = "hv")
  cowplot::plot_grid(title_grob, combined, ncol = 1, rel_heights = c(0.07, 1))
}


#' Save PPCA Plot for a Complex to PDF
#'
#' Wrapper around \code{plot_grouping_pca} that writes the combined plot to a
#' PDF file.
#'
#' @param vsclust_result Result from vsclust_on_complex.
#' @param intensity_cols Character vector of raw intensity column names.
#' @param n_rep Integer. Number of replicates per condition.
#' @param cond_regex Character. Regex for extracting condition labels.
#' @param complex_name Character. Complex name used for the plot title.
#' @param file_path Character. Output PDF path.
#' @param width Numeric. PDF width in inches (default 10).
#' @param height Numeric. PDF height in inches (default 5).
#'
#' @return Invisibly returns the combined plot object.
#'
#' @export
save_grouping_pca_pdf <- function(
  vsclust_result,
  intensity_cols,
  n_rep,
  cond_regex,
  complex_name = "",
  file_path,
  width = 10,
  height = 5
) {
  p <- plot_grouping_pca(
    vsclust_result = vsclust_result,
    intensity_cols = intensity_cols,
    n_rep = n_rep,
    cond_regex = cond_regex,
    complex_name = complex_name
  )

  grDevices::pdf(file = file_path, width = width, height = height, onefile = FALSE)
  print(p)
  grDevices::dev.off()

  invisible(p)
}
