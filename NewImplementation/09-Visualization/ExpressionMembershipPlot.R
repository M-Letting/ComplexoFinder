# Libraries
library(ComplexHeatmap)
library(ggplot2)
library(grid)
library(circlize)
library(data.table)

#' Create Combined Complexoform Visualization from Pipeline Results
#'
#' Wrapper function to create visualization from complex_obj$vsclust results.
#' The members of each cluster (membership >= membership_cutoff) are shown in
#' the cluster profiles, colored by membership.
#'
#' @param vsclust_result Result from vsclust_on_complex (complex_obj$vsclust[[complex]])
#' @param complex_name Character. Complex name for plot titles
#' @param membership_cutoff Numeric. Membership a feature needs to be shown as
#'   a member of its cluster (default: 0.5, as used for quantification)
#'
#' @return Invisibly returns list with plot objects
#'
#' @export
create_complexoform_plot_from_pipeline <- function(
  vsclust_result,
  complex_name,
  membership_cutoff = 0.5
) {
  create_combined_complexoform_plot(
    ClustOut = vsclust_result$ClustOut,
    title = complex_name,
    membership_cutoff = membership_cutoff
  )
}

#' Create Combined Complexoform Visualization
#'
#' Creates a combined visualization with cluster profile panels on the left
#' (one per cluster, lines colored by membership value) and dual heatmaps
#' (intensity + membership) on the right.
#'
#' Features with membership >= membership_cutoff are shown in cluster profiles
#' with color scale ranging from the cutoff to 1.
#'
#' @param ClustOut Clustering result from runClustWrapper_custom
#' @param title Character. plot titles
#' @param membership_cutoff Numeric. Membership a feature needs to be shown as
#'   a member of its cluster (default: 0.5)
#'
#' @return Invisibly returns list with plot objects
#'
#' @export
create_combined_complexoform_plot <- function(
  ClustOut,
  title,
  membership_cutoff = 0.5
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop(
      "Package 'ComplexHeatmap' is required. Install with: BiocManager::install('ComplexHeatmap')"
    )
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop(
      "Package 'circlize' is required. Install with: install.packages('circlize')"
    )
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop(
      "Package 'cowplot' is required. Install with: install.packages('cowplot')"
    )
  }

  # ---- Prepare data for profiles ----
  # One column per condition, named by its label (e.g., "D17", "D19", etc.)
  ht_data <- ClustOut$dat
  timepoint_names <- colnames(ht_data)

  n_clusters <- nrow(ClustOut$Bestcl$centers)
  cluster_assignments <- ClustOut$Bestcl$cluster
  membership_matrix <- ClustOut$Bestcl$membership

  # Ensure membership matrix is a matrix (not dropped to vector for k=1)
  if (!is.matrix(membership_matrix)) {
    membership_matrix <- matrix(membership_matrix, ncol = 1)
    rownames(membership_matrix) <- rownames(ht_data)
  }

  # Create long-format data for profile plots
  profile_long <- data.table(
    Feature = rep(rownames(ht_data), each = ncol(ht_data)),
    Timepoint = rep(1:ncol(ht_data), nrow(ht_data)),
    TimepointLabel = rep(timepoint_names, nrow(ht_data)),
    Value = as.vector(t(ht_data)),
    Cluster = rep(cluster_assignments, each = ncol(ht_data))
  )

  # Add membership values for each feature
  # Handle k=1 case where membership matrix has only 1 column
  profile_long[,
    Membership := {
      feature_idx <- match(Feature, rownames(ht_data))
      cluster_col <- min(Cluster[1], ncol(membership_matrix))
      membership_matrix[feature_idx, cluster_col]
    },
    by = Feature
  ]

  # ---- Create cluster profile plots ----
  cluster_plots <- list()

  # Calculate global y-axis limits across all clusters for consistent scaling
  y_min <- min(profile_long$Value, na.rm = TRUE)
  y_max <- max(profile_long$Value, na.rm = TRUE)

  for (k in 1:n_clusters) {
    cluster_data <- profile_long[Cluster == k]

    # Only show the members of the cluster
    cluster_data <- cluster_data[Membership >= membership_cutoff]

    n_features <- length(unique(cluster_data$Feature))

    p <- ggplot() +
      # Individual feature lines colored by membership (yellow -> orange -> red)
      geom_line(
        data = cluster_data,
        aes(x = Timepoint, y = Value, group = Feature, color = Membership),
        alpha = 0.7,
        linewidth = 0.5,
        na.rm = TRUE
      ) +
      scale_color_gradientn(
        colors = scico::scico(
          n = 5,
          palette = "lajolla",
          begin = 0.35,
          end = 1,
          direction = -1
        ),
        limits = c(membership_cutoff, 1),
        name = "Membership"
      ) +
      scale_x_continuous(
        breaks = seq(1, ncol(ht_data), by = max(1, floor(ncol(ht_data) / 7))),
        labels = timepoint_names[seq(
          1,
          ncol(ht_data),
          by = max(1, floor(ncol(ht_data) / 7))
        )]
      ) +
      coord_cartesian(ylim = c(y_min, y_max)) +
      labs(title = paste0("Cluster ", k), x = NULL, y = "Expression") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5)
      )

    cluster_plots[[k]] <- p
  }

  # Add a shared legend
  legend_plot <- ggplot(
    data.frame(x = seq(membership_cutoff, 1, length.out = 50)),
    aes(x, x, color = x)
  ) +
    geom_point() +
    scale_color_gradientn(
      colors = scico::scico(
        n = 5,
        palette = "lajolla",
        begin = 0.35,
        end = 1,
        direction = -1
      ),
      limits = c(membership_cutoff, 1),
      name = "Membership"
    ) +
    theme_void() +
    theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"))
  legend <- cowplot::get_legend(legend_plot)

  # Calculate number of columns for profile plots
  # 1-4 clusters: 1 column, 5-8: 2 columns, 9-12: 3 columns, etc.
  profile_ncol <- ceiling(n_clusters / 4)

  # Combine cluster plots with dynamic columns
  profiles_combined <- cowplot::plot_grid(
    plotlist = cluster_plots,
    ncol = profile_ncol,
    align = "hv"
  )
  profiles_with_legend <- cowplot::plot_grid(
    profiles_combined,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.05)
  )

  # ---- Create dual heatmaps ----
  # Split rownames by _
  rownames_split <- strsplit(rownames(ht_data), "_")

  ht_data_dt <- as.data.table(ht_data)
  ht_data_dt[, Gene_name := sapply(rownames_split, function(x) x[1])]
  ht_data_dt[, Accession := sapply(rownames_split, function(x) x[2])]
  ht_data_dt[,
    Position_in_master_protein := sapply(rownames_split, function(x) x[3])
  ]
  ht_data_dt[,
    Modifications_in_master_protein := sapply(rownames_split, function(x) x[4])
  ]
  ht_data_dt[, datatype := sapply(rownames_split, function(x) x[5])]
  ht_data_dt[, datatype := gsub("\\.\\d+$", "", datatype)]
  ht_data_dt[, Cluster := as.character(cluster_assignments)]

  max_membership <- sapply(1:nrow(ht_data_dt), function(i) {
    cluster_num <- as.numeric(ht_data_dt$Cluster[i])
    membership_matrix[i, cluster_num]
  })
  ht_data_dt[, max_membership := max_membership]

  # Order by cluster and membership
  row_order <- order(as.numeric(ht_data_dt$Cluster), -ht_data_dt$max_membership)
  ht_data_dt <- ht_data_dt[row_order, ]

  # Extract matrix for heatmap
  ht_matrix <- as.matrix(ht_data_dt[,
    !c(
      "Accession",
      "Gene_name",
      "Position_in_master_protein",
      "Modifications_in_master_protein",
      "datatype",
      "Cluster",
      "max_membership"
    ),
    with = FALSE
  ])
  ht_matrix[is.infinite(ht_matrix)] <- NA

  # Define colors
  datatype_colors <- c(
    "NM" = "#00441b",
    "NMpeptide" = "#02d054ff",
    "RmCys" = "#006d2c",
    "Deglyco" = "#99d38d",
    "Phospho" = "#be4e16",
    "FreeCys" = "#ff8000",
    "LysAc" = "#fec389"
  )

  cluster_colors <- scales::hue_pal()(n_clusters)
  names(cluster_colors) <- as.character(1:n_clusters)

  row_ha <- rowAnnotation(
    Datatype = ht_data_dt$datatype,
    col = list(Datatype = datatype_colors),
    show_legend = TRUE
  )

  intensity_range <- range(ht_matrix, na.rm = TRUE)
  intensity_breaks <- seq(
    intensity_range[1],
    intensity_range[2],
    length.out = 9
  )
  rownames(ht_matrix) <- ht_data_dt$Position_in_master_protein

  ht1 <- Heatmap(
    ht_matrix,
    name = "Intensity",
    col = colorRamp2(
      intensity_breaks,
      colorRampPalette(c(
        "#030b47",
        "#2778a3",
        "#e1e5e5",
        "#c47145",
        "#600002"
      ))(9)
    ),
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = NULL,
    right_annotation = row_ha
  )

  # Membership heatmap
  membership_matrix_ordered <- ClustOut$Bestcl$membership[
    row_order,
    ,
    drop = FALSE
  ]

  # Handle column naming - check if we have columns
  if (ncol(membership_matrix_ordered) > 0) {
    colnames(membership_matrix_ordered) <- paste0(
      "C",
      1:ncol(membership_matrix_ordered)
    )
  }

  membership_matrix_ordered[membership_matrix_ordered == 0] <- NA
  membership_range <- range(membership_matrix_ordered, na.rm = TRUE)

  # Handle case where all membership values are the same (e.g., single cluster)
  if (diff(membership_range) < 1e-10) {
    # All values are essentially the same - expand range slightly
    membership_breaks <- c(
      membership_range[1] - 0.1,
      membership_range[1],
      membership_range[1] + 0.1
    )
    membership_colors <- scico::scico(
      n = 3,
      palette = "lajolla",
      begin = 0.35,
      end = 1,
      direction = -1
    )
  } else {
    membership_breaks <- seq(
      membership_range[1],
      membership_range[2],
      length.out = 9
    )
    membership_colors <- scico::scico(
      n = length(membership_breaks),
      palette = "lajolla",
      begin = 0.35,
      end = 1,
      direction = -1
    )
  }

  row_ha_membership <- rowAnnotation(
    Datatype = ht_data_dt$datatype,
    col = list(Datatype = datatype_colors),
    show_legend = FALSE
  )

  ht2 <- Heatmap(
    membership_matrix_ordered,
    name = "Membership",
    col = colorRamp2(membership_breaks, membership_colors),
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = NULL,
    right_annotation = row_ha_membership
  )

  # ---- Combine: profiles on left, heatmaps on right ----
  # Draw into the currently active graphics device.

  # Create layout with title at top
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      2,
      1,
      heights = grid::unit(c(0.05, 0.95), "null")
    )
  ))

  # Draw main title
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.text(
    title,
    gp = grid::gpar(fontsize = 16, fontface = "bold"),
    y = 0.5
  )
  grid::popViewport()

  # Create horizontal layout for content: profiles (1/3) | heatmaps (2/3)
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(1, 2, widths = grid::unit(c(1, 2), "null"))
  ))

  # Draw profile plots in left panel
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(
    profiles_with_legend,
    vp = grid::viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.95)
  )
  grid::popViewport()

  # Draw heatmaps in right panel
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(ht1 + ht2, newpage = FALSE)
  grid::popViewport()

  grid::popViewport()
  grid::popViewport()

  invisible(list(cluster_plots = cluster_plots, ht1 = ht1, ht2 = ht2))
}
