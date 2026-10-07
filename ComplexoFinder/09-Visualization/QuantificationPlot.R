library(data.table)
library(ggplot2)
library(scico)

complexoform_palette <- function(n) {
  scico::scico(n, palette = "batlow")
}

# "membership of cluster N" -> "Complexoform group N"
clean_group_labels <- function(labels) {
  sub("^membership of cluster (\\d+)$", "Complexoform group \\1", labels)
}

#' Line plot of complexoform group abundance across conditions
#'
#' @param quantification_result List. Output of
#'   `quantify_complexoform_abundance_from_vsclust()` with `cond_regex` supplied
#'   so `abundance_wide` has one column per condition.
#' @param complex_name Character. Optional plot title.
#' @param group_colors Named character vector mapping group labels to colors.
#'   Auto-generated from `complexoform_palette()` when NULL.
#'
#' @return A ggplot object.
#'
#' @export
plot_complexoform_quantification <- function(
  quantification_result,
  complex_name = NULL,
  group_colors = NULL
) {
  abundance_wide <- copy(as.data.table(quantification_result$abundance_wide))
  condition_cols <- setdiff(names(abundance_wide), "ComplexoformGroup")

  long_dt <- data.table::melt(
    abundance_wide,
    id.vars = "ComplexoformGroup",
    measure.vars = condition_cols,
    variable.name = "Condition",
    value.name = "Abundance"
  )

  long_dt[, Group := clean_group_labels(ComplexoformGroup)]

  # Sort conditions numerically (D1, D2, ..., D10) with alphabetical fallback
  numeric_parts <- suppressWarnings(
    as.numeric(gsub("[^0-9]", "", condition_cols))
  )
  cond_order <- if (anyNA(numeric_parts)) {
    sort(condition_cols)
  } else {
    condition_cols[order(numeric_parts)]
  }
  long_dt[, Condition := factor(as.character(Condition), levels = cond_order)]

  groups <- unique(long_dt$Group)
  n_groups <- length(groups)

  if (is.null(group_colors)) {
    pal <- complexoform_palette(n_groups)
    names(pal) <- groups
  } else {
    pal <- group_colors
  }

  ggplot(long_dt, aes(x = Condition, y = Abundance, color = Group, group = Group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_color_manual(values = pal, name = "Complexoform group") +
    labs(title = complex_name, x = "Condition", y = "Abundance") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

#' Save complexoform quantification plot as PDF (16 x 6)
#'
#' @param quantification_result List. Output of
#'   `quantify_complexoform_abundance_from_vsclust()`.
#' @param file_path Character. Output PDF path.
#' @param complex_name Character. Optional plot title.
#'
#' @return Invisibly returns the plot object.
#'
#' @export
save_complexoform_quantification_plot <- function(
  quantification_result,
  file_path,
  complex_name = NULL
) {
  p <- plot_complexoform_quantification(
    quantification_result = quantification_result,
    complex_name = complex_name
  )

  grDevices::pdf(file = file_path, width = 16, height = 6, onefile = FALSE)
  print(p)
  grDevices::dev.off()

  invisible(p)
}
