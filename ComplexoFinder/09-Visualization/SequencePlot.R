library(data.table)
library(ggplot2)
library(cowplot)

#' Datatype color palette used in sequence views
#'
#' @return Named character vector of hex colors by datatype.
datatype_sequence_palette <- function() {
  c(
    "NM" = "#00441b",
    "NMpeptide" = "#02d054ff",
    "RmCys" = "#006d2c",
    "Deglyco" = "#99d38d",
    "Phospho" = "#be4e16",
    "FreeCys" = "#ff8000",
    "LysAc" = "#fec389"
  )
}

# Parse accession [start-end] annotations from a position string.
parse_position_annotations <- function(
  position_string,
  accession_hint = NA_character_
) {
  if (is.na(position_string) || !nzchar(trimws(position_string))) {
    return(data.table(
      Accession = character(0),
      start_position = integer(0),
      end_position = integer(0)
    ))
  }

  tokens <- trimws(strsplit(position_string, ";", fixed = TRUE)[[1]])
  current_accession <- accession_hint
  parsed <- vector("list", length(tokens))
  out_idx <- 0L

  for (tok in tokens) {
    if (!nzchar(tok)) {
      next
    }

    if (grepl("^[^\\[]+\\s*\\[", tok)) {
      current_accession <- trimws(sub("\\s*\\[.*$", "", tok))
    }

    if (grepl("\\[NOT FOUND\\]", tok)) {
      next
    }

    ranges <- regmatches(tok, gregexpr("\\[\\d+-\\d+\\]", tok))[[1]]
    if (length(ranges) == 0) {
      next
    }

    for (rg in ranges) {
      start_val <- as.integer(sub("\\[(\\d+)-\\d+\\]", "\\1", rg))
      end_val <- as.integer(sub("\\[\\d+-(\\d+)\\]", "\\1", rg))
      if (is.na(start_val) || is.na(end_val)) {
        next
      }

      out_idx <- out_idx + 1L
      parsed[[out_idx]] <- data.table(
        Accession = as.character(current_accession),
        start_position = start_val,
        end_position = end_val
      )
    }
  }

  if (out_idx == 0L) {
    return(data.table(
      Accession = character(0),
      start_position = integer(0),
      end_position = integer(0)
    ))
  }

  rbindlist(parsed[seq_len(out_idx)])
}

# Compare accessions with isoform-1 normalization.
# Isoform-1 (X-1) is treated as equivalent to the canonical base (X).
# Isoforms 2+ require an exact match, so Q13936-2 != Q13936.
same_accession <- function(x, y) {
  norm <- function(a) sub("-1$", "", as.character(a))
  xn <- norm(x)
  yn <- norm(y)
  !is.na(xn) & !is.na(yn) & xn == yn
}

# Fetch protein sequence using existing get_sequence() when available.
fetch_sequence_for_accession <- function(accession) {
  accession <- as.character(accession)
  candidates <- unique(c(accession, sub("-\\d+$", "", accession)))

  for (acc in candidates) {
    seq_try <- tryCatch(
      {
        if (exists("get_sequence", mode = "function")) {
          get_sequence(acc)$sequence
        } else {
          fasta <- readLines(
            paste0("https://rest.uniprot.org/uniprotkb/", acc, ".fasta"),
            warn = FALSE
          )
          fasta_text <- paste(fasta, collapse = "\n")
          seq <- gsub(">[^\\n]+", "", fasta_text)
          gsub("\\n", "", seq)
        }
      },
      error = function(e) {
        NA_character_
      }
    )

    if (!is.na(seq_try) && nzchar(seq_try)) {
      return(seq_try)
    }
  }

  NA_character_
}

# Stack peptide intervals greedily to avoid overlap.
stack_intervals <- function(start_pos, end_pos) {
  if (length(start_pos) == 0L) {
    return(integer(0))
  }

  last_end <- numeric(0)
  level_idx <- integer(length(start_pos))

  for (i in seq_along(start_pos)) {
    s <- start_pos[i]
    e <- end_pos[i]
    candidate <- which(s > last_end)

    if (length(candidate) > 0) {
      j <- candidate[1]
      last_end[j] <- max(last_end[j], e)
    } else {
      j <- length(last_end) + 1L
      last_end[j] <- e
    }

    level_idx[i] <- j
  }

  level_idx
}

# Build one sequence panel for a single accession.
create_accession_datatype_sequence_plot <- function(
  accession_data,
  accession,
  sequence,
  complex_name = NULL,
  datatype_col = "datatype"
) {
  y_height <- 1
  y_offset <- 1.1
  y_step <- 1.1

  seq_len <- if (!is.na(sequence) && nzchar(sequence)) {
    nchar(sequence)
  } else {
    max(accession_data$end_position, na.rm = TRUE)
  }

  x_step <- max(50, ceiling(seq_len / 10 / 50) * 50)
  palette <- datatype_sequence_palette()

  accession_data <- copy(accession_data)
  accession_data <- accession_data[order(start_position, end_position)]

  is_nmpeptide <- accession_data[[datatype_col]] == "NMpeptide"
  ptm_data <- accession_data[!is_nmpeptide]
  nmpeptide_data <- accession_data[is_nmpeptide]

  if (nrow(ptm_data) > 0) {
    ptm_levels <- stack_intervals(
      ptm_data$start_position,
      ptm_data$end_position
    )
    ptm_data[, y_center := y_offset + (ptm_levels - 1L) * y_step]
  }

  if (nrow(nmpeptide_data) > 0) {
    nm_levels <- stack_intervals(
      nmpeptide_data$start_position,
      nmpeptide_data$end_position
    )
    nmpeptide_data[, y_center := -y_offset - (nm_levels - 1L) * y_step]
  }

  rect_data <- rbindlist(
    list(ptm_data, nmpeptide_data),
    use.names = TRUE,
    fill = TRUE
  )
  if (nrow(rect_data) > 0) {
    rect_data[, ymin := y_center - y_height / 2]
    rect_data[, ymax := y_center + y_height / 2]
  }

  plot_title <- if (!is.null(complex_name) && nzchar(complex_name)) {
    paste0(complex_name, " | ", accession)
  } else {
    accession
  }

  p <- ggplot() +
    geom_rect(
      data = data.frame(xmin = 1, xmax = seq_len, y_level = 0),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = y_level - y_height / 2,
        ymax = y_level + y_height / 2
      ),
      fill = palette[["NM"]],
      alpha = 1
    ) +
    labs(title = plot_title, x = "Amino acid position", y = NULL) +
    scale_x_continuous(
      breaks = seq(0, seq_len, by = x_step),
      expand = c(0.01, 0.01)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "none"
    )

  if (nrow(rect_data) > 0) {
    p <- p +
      geom_rect(
        data = rect_data,
        aes(
          xmin = start_position,
          xmax = end_position,
          ymin = ymin,
          ymax = ymax,
          fill = .data[[datatype_col]]
        ),
        color = "black",
        linewidth = 0.25,
        alpha = 1
      ) +
      scale_fill_manual(values = palette, drop = FALSE)

    max_y <- max(rect_data$ymax, na.rm = TRUE)
    min_y <- min(rect_data$ymin, na.rm = TRUE)
    p <- p + coord_cartesian(ylim = c(min(min_y, -1), max(max_y, 1)))
  }

  p
}

#' Plot datatype-aware sequence views for each accession in a complex
#'
#' Creates one sequence plot per accession found in complex peptide data and
#' combines all panels into a single cowplot object with a shared datatype legend.
#'
#' @param complex_data data.frame/data.table with at least `Accession`, datatype,
#'   and peptide position annotations.
#' @param complex_name Character. Optional title prefix for each panel.
#' @param accession_col Character. Accession column name.
#' @param datatype_col Character. Datatype column name.
#' @param position_col Character. Position annotation column. If `NULL`, auto-detects
#'   `"Position in master protein"` or `"Position"`.
#' @param ncol Integer. Number of columns used for accession panels.
#'
#' @return A combined `ggplot` object created with `cowplot::plot_grid`.
#'
#' @export
plot_sequence_views_datatype_pipeline <- function(
  complex_data,
  complex_name = NULL,
  accession_col = "Accession",
  datatype_col = "datatype",
  position_col = NULL,
  ncol = 1
) {
  if (!is.data.table(complex_data)) {
    complex_data <- as.data.table(complex_data)
  }

  required_cols <- c(accession_col, datatype_col)
  missing_required <- setdiff(required_cols, names(complex_data))
  if (length(missing_required) > 0) {
    stop("Missing required columns: ", paste(missing_required, collapse = ", "))
  }

  if (is.null(position_col)) {
    if ("Position in master protein" %in% names(complex_data)) {
      position_col <- "Position in master protein"
    } else if ("Position" %in% names(complex_data)) {
      position_col <- "Position"
    } else {
      stop("No position column found. Provide position_col explicitly.")
    }
  }

  raw_accessions <- unique(unlist(strsplit(
    complex_data[[accession_col]],
    ";",
    fixed = TRUE
  )))
  accessions <- trimws(raw_accessions)
  accessions <- accessions[!is.na(accessions) & nzchar(accessions)]
  accessions <- unique(accessions)

  if (length(accessions) == 0) {
    stop("No accession values found in complex_data.")
  }

  panel_plots <- list()
  all_dtypes <- character(0)

  for (acc in accessions) {
    mapped_rows <- vector("list", nrow(complex_data))
    mapped_idx <- 0L

    for (i in seq_len(nrow(complex_data))) {
      row_acc <- as.character(complex_data[[accession_col]][i])
      hint <- trimws(strsplit(row_acc, ";", fixed = TRUE)[[1]])
      hint <- hint[!is.na(hint) & nzchar(hint)]
      first_hint <- if (length(hint) > 0) hint[1] else NA_character_

      parsed <- parse_position_annotations(
        position_string = as.character(complex_data[[position_col]][i]),
        accession_hint = first_hint
      )

      if (nrow(parsed) == 0) {
        next
      }

      keep <- same_accession(parsed$Accession, acc)

      if (!any(keep) && length(hint) == 1 && same_accession(hint, acc)) {
        keep <- rep(TRUE, nrow(parsed))
        parsed[, Accession := hint]
      }

      if (!any(keep)) {
        next
      }

      parsed <- parsed[keep]
      parsed[, (datatype_col) := complex_data[[datatype_col]][i]]

      mapped_idx <- mapped_idx + 1L
      mapped_rows[[mapped_idx]] <- parsed
    }

    if (mapped_idx == 0L) {
      next
    }

    accession_data <- rbindlist(
      mapped_rows[seq_len(mapped_idx)],
      use.names = TRUE,
      fill = TRUE
    )
    accession_data <- accession_data[
      !is.na(start_position) & !is.na(end_position)
    ]

    if (nrow(accession_data) == 0) {
      next
    }

    sequence <- fetch_sequence_for_accession(acc)
    panel_plots[[acc]] <- create_accession_datatype_sequence_plot(
      accession_data = accession_data,
      accession = acc,
      sequence = sequence,
      complex_name = complex_name,
      datatype_col = datatype_col
    )

    all_dtypes <- union(
      all_dtypes,
      unique(as.character(accession_data[[datatype_col]]))
    )
  }

  if (length(panel_plots) == 0) {
    stop("No accession-specific sequence panels could be created.")
  }

  palette <- datatype_sequence_palette()
  legend_types <- setdiff(names(palette)[names(palette) %in% all_dtypes], "NM")

  main_panels <- cowplot::plot_grid(
    plotlist = panel_plots,
    ncol = ncol,
    align = "v"
  )

  if (length(legend_types) > 0) {
    legend_df <- data.frame(
      datatype = factor(legend_types, levels = legend_types),
      x = 1,
      y = 1
    )

    legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = datatype)) +
      geom_tile() +
      scale_fill_manual(
        values = palette[legend_types],
        breaks = legend_types,
        name = "Datatype",
        drop = FALSE
      ) +
      theme_void() +
      theme(legend.position = "right")

    shared_legend <- cowplot::get_legend(legend_plot)

    combined_plot <- cowplot::plot_grid(
      main_panels,
      shared_legend,
      ncol = 2,
      rel_widths = c(9, 1)
    )
  } else {
    combined_plot <- main_panels
  }

  attr(combined_plot, "n_accessions") <- length(panel_plots)
  attr(combined_plot, "legend_datatypes") <- legend_types

  combined_plot
}

# ---- Cluster membership sequence plots ----------------------------------------

# Build one sequence panel per accession, colored by cluster membership.
create_accession_cluster_sequence_plot <- function(
  accession_data,
  accession,
  sequence,
  cluster_label,
  datatype_col = "datatype",
  member_color = "#78185E",
  nonmember_color = "#969696"
) {
  y_height <- 1
  y_offset <- 1.1
  y_step <- 1.1

  seq_len_val <- if (!is.na(sequence) && nzchar(sequence)) {
    nchar(sequence)
  } else {
    max(accession_data$end_position, na.rm = TRUE)
  }

  x_step <- max(50, ceiling(seq_len_val / 10 / 50) * 50)

  accession_data <- copy(accession_data)
  accession_data[,
    membership_status := ifelse(is_member, "Member", "Non-member")
  ]
  accession_data <- accession_data[order(start_position, end_position)]

  is_nmpeptide <- accession_data[[datatype_col]] == "NMpeptide"
  ptm_data <- accession_data[!is_nmpeptide]
  nmpeptide_data <- accession_data[is_nmpeptide]

  if (nrow(ptm_data) > 0) {
    ptm_levels <- stack_intervals(
      ptm_data$start_position,
      ptm_data$end_position
    )
    ptm_data[, y_center := y_offset + (ptm_levels - 1L) * y_step]
  }
  if (nrow(nmpeptide_data) > 0) {
    nm_levels <- stack_intervals(
      nmpeptide_data$start_position,
      nmpeptide_data$end_position
    )
    nmpeptide_data[, y_center := -y_offset - (nm_levels - 1L) * y_step]
  }

  rect_data <- rbindlist(
    list(ptm_data, nmpeptide_data),
    use.names = TRUE,
    fill = TRUE
  )
  if (nrow(rect_data) > 0) {
    rect_data[, ymin := y_center - y_height / 2]
    rect_data[, ymax := y_center + y_height / 2]
  }

  p <- ggplot() +
    geom_rect(
      data = data.frame(xmin = 1, xmax = seq_len_val, y_level = 0),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = y_level - y_height / 2,
        ymax = y_level + y_height / 2
      ),
      fill = "#00441b",
      alpha = 1
    ) +
    labs(
      title = paste0(cluster_label, " | ", accession),
      x = "Amino acid position",
      y = NULL
    ) +
    scale_x_continuous(
      breaks = seq(0, seq_len_val, by = x_step),
      expand = c(0.01, 0.01)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "none"
    )

  if (nrow(rect_data) > 0) {
    p <- p +
      geom_rect(
        data = rect_data,
        aes(
          xmin = start_position,
          xmax = end_position,
          ymin = ymin,
          ymax = ymax,
          fill = membership_status
        ),
        color = "black",
        linewidth = 0.25,
        alpha = 1
      ) +
      scale_fill_manual(
        values = c("Member" = member_color, "Non-member" = nonmember_color),
        drop = FALSE
      )

    max_y <- max(rect_data$ymax, na.rm = TRUE)
    min_y <- min(rect_data$ymin, na.rm = TRUE)
    p <- p + coord_cartesian(ylim = c(min(min_y, -1), max(max_y, 1)))
  }

  p
}

#' Plot cluster-membership sequence views for each accession in a complex
#'
#' Creates one combined sequence plot per cluster, with one panel per accession.
#' Peptides with membership >= 1/k (k = number of clusters) are colored in
#' member_color; all others use nonmember_color. Title for each cluster plot is
#' "Complexoform group n".
#'
#' @param vsclust_result List. Output of `vsclust_on_complex()`.
#' @param complex_data data.frame/data.table with Accession, datatype, and
#'   position annotation columns (the original per-complex data, not the
#'   stripped peptide_results table).
#' @param id_cols Character vector. Column names used to construct VSClust row
#'   identifiers; must match those passed to `vsclust_on_complex()`.
#' @param accession_col Character. Accession column name in `complex_data`.
#' @param datatype_col Character. Datatype column name in `complex_data`.
#' @param position_col Character. Position annotation column. Auto-detected if NULL.
#' @param ncol Integer. Number of panel columns per cluster figure.
#' @param member_color Character. Hex color for member peptides.
#' @param nonmember_color Character. Hex color for non-member peptides.
#'
#' @return Named list of combined cowplot objects (one per cluster, named
#'   "cluster_1", "cluster_2", ...). Each carries an `n_accessions` attribute
#'   for dynamic height calculation.
#'
#' @export
plot_sequence_views_cluster_pipeline <- function(
  vsclust_result,
  complex_data,
  id_cols = c("Gene name", "Peptide"),
  accession_col = "Accession",
  datatype_col = "datatype",
  position_col = NULL,
  ncol = 1,
  member_color = "#78185E",
  nonmember_color = "#969696"
) {
  if (!is.data.table(complex_data)) {
    complex_data <- as.data.table(complex_data)
  }

  if (is.null(position_col)) {
    if ("Position in master protein" %in% names(complex_data)) {
      position_col <- "Position in master protein"
    } else if ("Position" %in% names(complex_data)) {
      position_col <- "Position"
    } else {
      stop("No position column found. Provide position_col explicitly.")
    }
  }

  membership_mat <- vsclust_result$ClustOut$Bestcl$membership
  k <- ncol(membership_mat)
  threshold <- 1 / 2

  # Build lookup keys for complex_data rows to match membership matrix row names.
  # complex_data may lack a synthetic "Peptide" column that discover_complexoforms
  # creates from Accession + Position + Modifications + datatype. Reconstruct it
  # from those source columns when needed. make.unique suffixes (.1, .2, ...) are
  # stripped later from membership row names so duplicates are handled by max().
  if (all(id_cols %in% names(complex_data))) {
    data_keys <- apply(complex_data[, id_cols, with = FALSE], 1, function(row) {
      paste(row, collapse = "_")
    })
  } else {
    composite_src_cols <- c(
      "Accession",
      "Position in master protein",
      "Modifications in master protein",
      "datatype"
    )
    if (!all(composite_src_cols %in% names(complex_data))) {
      stop(
        "Cannot map complex_data rows to membership matrix. ",
        "Ensure complex_data contains id_cols (",
        paste(id_cols, collapse = ", "),
        ") or the composite source columns: ",
        paste(composite_src_cols, collapse = ", ")
      )
    }
    gene_part <- as.character(complex_data[["Gene name"]])
    peptide_part <- apply(
      complex_data[, composite_src_cols, with = FALSE],
      1,
      function(row) paste(row, collapse = "_")
    )
    data_keys <- paste(gene_part, peptide_part, sep = "_")
  }

  mem_rownames <- rownames(membership_mat)
  mem_base_keys <- sub("\\.[0-9]+$", "", mem_rownames)

  raw_accessions <- unique(unlist(strsplit(
    complex_data[[accession_col]],
    ";",
    fixed = TRUE
  )))
  accessions <- unique(trimws(raw_accessions[
    !is.na(raw_accessions) & nzchar(trimws(raw_accessions))
  ]))

  if (length(accessions) == 0) {
    stop("No accession values found in complex_data.")
  }

  sequences <- setNames(
    lapply(accessions, fetch_sequence_for_accession),
    accessions
  )

  cluster_plots <- vector("list", k)
  names(cluster_plots) <- paste0("cluster_", seq_len(k))

  for (cl in seq_len(k)) {
    cluster_label <- paste("Complexoform group", cl)
    mem_col <- paste("membership of cluster", cl)

    # For each complex_data row, take the max membership across any matching
    # vsclust rows (handles duplicates from make.unique).
    is_member_vec <- vapply(
      data_keys,
      function(key) {
        matches <- which(mem_base_keys == key)
        if (length(matches) == 0L) {
          return(FALSE)
        }
        vals <- membership_mat[matches, mem_col, drop = TRUE]
        !all(is.na(vals)) && max(vals, na.rm = TRUE) >= threshold
      },
      logical(1)
    )

    panel_plots <- list()

    for (acc in accessions) {
      mapped_rows <- vector("list", nrow(complex_data))
      mapped_idx <- 0L

      for (i in seq_len(nrow(complex_data))) {
        row_acc <- as.character(complex_data[[accession_col]][i])
        hint <- trimws(strsplit(row_acc, ";", fixed = TRUE)[[1]])
        hint <- hint[!is.na(hint) & nzchar(hint)]
        first_hint <- if (length(hint) > 0) hint[1] else NA_character_

        parsed <- parse_position_annotations(
          position_string = as.character(complex_data[[position_col]][i]),
          accession_hint = first_hint
        )
        if (nrow(parsed) == 0) {
          next
        }

        keep <- same_accession(parsed$Accession, acc)
        if (!any(keep) && length(hint) == 1 && same_accession(hint, acc)) {
          keep <- rep(TRUE, nrow(parsed))
          parsed[, Accession := hint]
        }
        if (!any(keep)) {
          next
        }

        parsed <- parsed[keep]
        parsed[, (datatype_col) := complex_data[[datatype_col]][i]]
        parsed[, is_member := is_member_vec[i]]

        mapped_idx <- mapped_idx + 1L
        mapped_rows[[mapped_idx]] <- parsed
      }

      if (mapped_idx == 0L) {
        next
      }

      accession_data <- rbindlist(
        mapped_rows[seq_len(mapped_idx)],
        use.names = TRUE,
        fill = TRUE
      )
      accession_data <- accession_data[
        !is.na(start_position) & !is.na(end_position)
      ]
      if (nrow(accession_data) == 0) {
        next
      }

      panel_plots[[acc]] <- create_accession_cluster_sequence_plot(
        accession_data = accession_data,
        accession = acc,
        sequence = sequences[[acc]],
        cluster_label = cluster_label,
        datatype_col = datatype_col,
        member_color = member_color,
        nonmember_color = nonmember_color
      )
    }

    if (length(panel_plots) == 0) {
      warning(paste("No panels created for cluster", cl))
      next
    }

    main_panels <- cowplot::plot_grid(
      plotlist = panel_plots,
      ncol = ncol,
      align = "v"
    )

    legend_df <- data.frame(
      Status = factor(
        c("Member", "Non-member"),
        levels = c("Member", "Non-member")
      ),
      x = 1,
      y = 1
    )
    legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = Status)) +
      geom_tile() +
      scale_fill_manual(
        values = c("Member" = member_color, "Non-member" = nonmember_color),
        name = "Membership",
        drop = FALSE
      ) +
      theme_void() +
      theme(legend.position = "right")
    shared_legend <- cowplot::get_legend(legend_plot)

    combined_plot <- cowplot::plot_grid(
      main_panels,
      shared_legend,
      ncol = 2,
      rel_widths = c(9, 1)
    )
    attr(combined_plot, "n_accessions") <- length(panel_plots)
    cluster_plots[[paste0("cluster_", cl)]] <- combined_plot
  }

  cluster_plots[!vapply(cluster_plots, is.null, logical(1))]
}

#' Save cluster membership sequence plots as PDFs
#'
#' Saves one PDF per cluster to `output_dir`, named
#' `cluster_sequence_plot_<n>.pdf`. Height is dynamic: 2 * number of accessions.
#'
#' @param vsclust_result List. Output of `vsclust_on_complex()`.
#' @param complex_data data.frame/data.table with annotation columns (Accession,
#'   datatype, position).
#' @param output_dir Character. Directory to write PDF files into.
#' @param id_cols Character vector. Column names used to construct VSClust row
#'   identifiers.
#' @param accession_col Character. Accession column name.
#' @param datatype_col Character. Datatype column name.
#' @param position_col Character. Position annotation column.
#' @param ncol Integer. Number of panel columns per cluster figure.
#' @param member_color Character. Hex color for member peptides.
#' @param nonmember_color Character. Hex color for non-member peptides.
#'
#' @return Invisibly returns the named list of cluster plots.
#'
#' @export
save_cluster_sequence_plots_pdf <- function(
  vsclust_result,
  complex_data,
  output_dir,
  id_cols = c("Gene name", "Peptide"),
  accession_col = "Accession",
  datatype_col = "datatype",
  position_col = NULL,
  ncol = 1,
  member_color = "#78185E",
  nonmember_color = "#969696"
) {
  cluster_plots <- plot_sequence_views_cluster_pipeline(
    vsclust_result = vsclust_result,
    complex_data = complex_data,
    id_cols = id_cols,
    accession_col = accession_col,
    datatype_col = datatype_col,
    position_col = position_col,
    ncol = ncol,
    member_color = member_color,
    nonmember_color = nonmember_color
  )

  for (nm in names(cluster_plots)) {
    cl_num <- sub("^cluster_", "", nm)
    plot_obj <- cluster_plots[[nm]]

    n_acc <- attr(plot_obj, "n_accessions")
    if (is.null(n_acc) || is.na(n_acc) || n_acc < 1) {
      n_acc <- 1
    }

    grDevices::pdf(
      file = file.path(
        output_dir,
        paste0("cluster_sequence_plot_", cl_num, ".pdf")
      ),
      width = 16,
      height = 2 * n_acc,
      onefile = FALSE
    )
    print(plot_obj)
    grDevices::dev.off()
  }

  invisible(cluster_plots)
}

#' Save datatype sequence views for a complex as PDF
#'
#' Saves the combined accession panels to a PDF with fixed width 16 and
#' dynamic height equal to `2 * number_of_accessions`.
#'
#' @param complex_data data.frame/data.table with complex peptide annotations.
#' @param file_path Character. Output PDF path.
#' @param complex_name Character. Optional complex name for panel titles.
#' @param ncol Integer. Number of panel columns in the saved figure.
#'
#' @return Invisibly returns the combined plot.
#'
#' @export
save_datatype_sequence_plot_pdf <- function(
  complex_data,
  file_path,
  complex_name = NULL,
  ncol = 1
) {
  plot_obj <- plot_sequence_views_datatype_pipeline(
    complex_data = complex_data,
    complex_name = complex_name,
    ncol = ncol
  )

  n_accessions <- attr(plot_obj, "n_accessions")
  if (is.null(n_accessions) || is.na(n_accessions) || n_accessions < 1) {
    n_accessions <- 1
  }

  grDevices::pdf(
    file = file_path,
    width = 16,
    height = 2 * n_accessions,
    onefile = FALSE
  )
  print(plot_obj)
  grDevices::dev.off()

  invisible(plot_obj)
}
