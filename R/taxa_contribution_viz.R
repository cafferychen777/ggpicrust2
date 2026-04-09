# =============================================================================
# Taxa-Function Contribution Visualization
# =============================================================================
# Stacked bar plots and heatmaps showing which taxa drive pathway abundances.

#' Stacked bar plot of taxa contributions
#'
#' Creates a stacked bar plot showing taxa contributions to predicted
#' functional abundances, faceted by function or sample group.
#'
#' @param contrib_agg A data.frame from \code{\link{aggregate_taxa_contributions}}.
#' @param metadata A data.frame containing sample metadata.
#' @param group Character. Column name in \code{metadata} for grouping samples.
#' @param function_ids Optional character vector of function IDs to plot.
#'   If NULL (default), the top \code{n_functions} by variance are shown.
#' @param n_functions Integer. Number of functions to show when
#'   \code{function_ids} is NULL. Default 6.
#' @param facet_by Character. Facet by \code{"function"} (default) or
#'   \code{"group"}.
#' @param show_percentage Logical. Normalize bars to 100\%? Default TRUE.
#' @param color_theme Character. Color theme name, passed to
#'   \code{\link{get_color_theme}}. Default \code{"default"}.
#' @param font_size Numeric. Base font size. Default 12.
#' @param legend_position Character. Legend position. Default \code{"right"}.
#' @param custom_title Optional character string for the plot title.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \donttest{
#' # Synthetic example
#' agg <- data.frame(
#'   sample = rep(c("S1", "S2", "S3", "S4"), each = 3),
#'   function_id = rep(c("K00001", "K00002"), each = 6),
#'   taxon_label = rep(c("Genus_A", "Genus_B", "Other"), 4),
#'   contribution = runif(12)
#' )
#' metadata <- data.frame(
#'   sample = c("S1", "S2", "S3", "S4"),
#'   group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' p <- taxa_contribution_bar(agg, metadata, group = "group")
#' }
#'
#' @export
taxa_contribution_bar <- function(contrib_agg,
                                  metadata,
                                  group,
                                  function_ids = NULL,
                                  n_functions = 6,
                                  facet_by = "function",
                                  show_percentage = TRUE,
                                  color_theme = "default",
                                  font_size = 12,
                                  legend_position = "right",
                                  custom_title = NULL) {
  # Validate inputs
  validate_dataframe(contrib_agg,
                     required_cols = c("sample", "function_id",
                                       "taxon_label", "contribution"),
                     param_name = "contrib_agg")
  validate_metadata(metadata)
  validate_group(metadata, group, min_groups = 1)

  facet_by <- match.arg(facet_by, c("function", "group"))

  # Align samples
  # Build a pseudo-abundance matrix for align_samples()
  sample_ids <- unique(contrib_agg$sample)
  pseudo_abundance <- matrix(1, nrow = 1, ncol = length(sample_ids))
  colnames(pseudo_abundance) <- sample_ids
  rownames(pseudo_abundance) <- "dummy"
  aligned <- align_samples(pseudo_abundance, metadata, verbose = FALSE)
  common_samples <- colnames(aligned$abundance)
  metadata_aligned <- aligned$metadata

  contrib_agg <- contrib_agg[contrib_agg$sample %in% common_samples, ]
  if (nrow(contrib_agg) == 0) {
    stop("No overlapping samples between contrib_agg and metadata.")
  }

  # Select functions to plot
  if (is.null(function_ids)) {
    func_var <- stats::aggregate(
      contribution ~ function_id, data = contrib_agg, FUN = stats::var
    )
    func_var <- func_var[order(-func_var$contribution), ]
    function_ids <- utils::head(func_var$function_id, n_functions)
  }
  contrib_agg <- contrib_agg[contrib_agg$function_id %in% function_ids, ]

  # Add group info
  sample_col <- aligned$sample_col
  if (sample_col == ".rownames") {
    contrib_agg$group_var <- metadata_aligned[contrib_agg$sample, group]
  } else {
    group_map <- stats::setNames(
      metadata_aligned[[group]],
      metadata_aligned[[sample_col]]
    )
    contrib_agg$group_var <- group_map[contrib_agg$sample]
  }

  # Normalize to percentage if requested
  if (show_percentage) {
    contrib_agg <- do.call(rbind, lapply(
      split(contrib_agg, list(contrib_agg$sample, contrib_agg$function_id)),
      function(df) {
        total <- sum(df$contribution, na.rm = TRUE)
        if (total > 0) df$contribution <- df$contribution / total * 100
        df
      }
    ))
    rownames(contrib_agg) <- NULL
    y_label <- "Relative contribution (%)"
  } else {
    y_label <- "Contribution"
  }

  # Get colors
  n_taxa <- length(unique(contrib_agg$taxon_label))
  theme_colors <- get_color_theme(color_theme, n_colors = n_taxa)
  taxa_colors <- theme_colors$pathway_class_colors[seq_len(n_taxa)]
  names(taxa_colors) <- sort(unique(contrib_agg$taxon_label))
  # Ensure "Other" is grey
  if ("Other" %in% names(taxa_colors)) {
    taxa_colors[["Other"]] <- "#999999"
  }

  # Build plot
  p <- ggplot2::ggplot(
    contrib_agg,
    ggplot2::aes(x = .data$sample, y = .data$contribution,
                 fill = .data$taxon_label)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = taxa_colors, name = "Taxon") +
    ggplot2::labs(
      x = "Sample",
      y = y_label,
      title = custom_title
    ) +
    ggprism::theme_prism(base_size = font_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = legend_position
    )

  # Faceting
  if (facet_by == "function") {
    p <- p + ggplot2::facet_wrap(~ function_id, scales = "free_x")
  } else {
    p <- p + ggplot2::facet_wrap(~ group_var, scales = "free_x")
  }

  p
}


#' Heatmap of taxa contributions across functions
#'
#' Creates a heatmap showing mean taxa contributions across pathways/functions,
#' with optional clustering and pathway annotations.
#'
#' @param contrib_agg A data.frame from \code{\link{aggregate_taxa_contributions}}.
#' @param annotation_data Optional data.frame from \code{\link{pathway_annotation}}
#'   for replacing function IDs with readable descriptions.
#' @param n_functions Integer. Number of functions to include. Default 20.
#' @param cluster_rows Logical. Cluster rows (taxa)? Default TRUE.
#' @param cluster_cols Logical. Cluster columns (functions)? Default TRUE.
#' @param clustering_method Character. Method for \code{hclust}. Default \code{"complete"}.
#' @param clustering_distance Character. Distance metric. Default \code{"euclidean"}.
#' @param low_color Character. Color for low values. Default \code{"#f7f7f7"}.
#' @param high_color Character. Color for high values. Default \code{"#ca0020"}.
#' @param font_size Numeric. Base font size. Default 12.
#' @param dendro_line_size Numeric. Dendrogram line width. Default 0.5.
#' @param custom_title Optional plot title.
#'
#' @return A \code{ggplot2} or \code{patchwork} object.
#'
#' @examples
#' \donttest{
#' agg <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 6),
#'   function_id = rep(rep(c("K00001", "K00002", "K00003"), each = 2), 2),
#'   taxon_label = rep(c("Genus_A", "Genus_B"), 6),
#'   contribution = runif(12)
#' )
#' p <- taxa_contribution_heatmap(agg)
#' }
#'
#' @export
taxa_contribution_heatmap <- function(contrib_agg,
                                      annotation_data = NULL,
                                      n_functions = 20,
                                      cluster_rows = TRUE,
                                      cluster_cols = TRUE,
                                      clustering_method = "complete",
                                      clustering_distance = "euclidean",
                                      low_color = "#f7f7f7",
                                      high_color = "#ca0020",
                                      font_size = 12,
                                      dendro_line_size = 0.5,
                                      custom_title = NULL) {
  validate_dataframe(contrib_agg,
                     required_cols = c("sample", "function_id",
                                       "taxon_label", "contribution"),
                     param_name = "contrib_agg")

  # Compute mean contribution per taxon x function
  mean_contrib <- stats::aggregate(
    contribution ~ function_id + taxon_label,
    data = contrib_agg,
    FUN = mean
  )

  # Select top functions by total contribution
  func_totals <- stats::aggregate(
    contribution ~ function_id, data = mean_contrib, FUN = sum
  )
  func_totals <- func_totals[order(-func_totals$contribution), ]
  top_funcs <- utils::head(func_totals$function_id, n_functions)
  mean_contrib <- mean_contrib[mean_contrib$function_id %in% top_funcs, ]

  # Pivot to wide matrix: taxon_label (rows) x function_id (cols)
  mat <- tidyr::pivot_wider(
    mean_contrib,
    names_from = "function_id",
    values_from = "contribution",
    values_fill = 0
  )
  taxa_names <- mat$taxon_label
  mat <- as.matrix(mat[, -1, drop = FALSE])
  rownames(mat) <- taxa_names

  # Replace function IDs with annotations if available
  if (!is.null(annotation_data)) {
    desc_map <- stats::setNames(
      annotation_data$description,
      annotation_data$pathway
    )
    new_names <- desc_map[colnames(mat)]
    # Only replace where we found a match, truncate long names
    found <- !is.na(new_names)
    new_names[found] <- substr(new_names[found], 1, 50)
    new_names[!found] <- colnames(mat)[!found]
    colnames(mat) <- new_names
  }

  # Clustering
  row_order <- seq_len(nrow(mat))
  col_order <- seq_len(ncol(mat))
  row_hc <- NULL
  col_hc <- NULL

  if (cluster_rows && nrow(mat) > 1) {
    row_hc <- stats::hclust(stats::dist(mat, method = clustering_distance),
                            method = clustering_method)
    row_order <- row_hc$order
  }
  if (cluster_cols && ncol(mat) > 1) {
    col_hc <- stats::hclust(stats::dist(t(mat), method = clustering_distance),
                            method = clustering_method)
    col_order <- col_hc$order
  }

  # Reorder matrix
  mat <- mat[row_order, col_order, drop = FALSE]

  # Convert to long for ggplot
  plot_df <- data.frame(
    taxon = rep(rownames(mat), ncol(mat)),
    func = rep(colnames(mat), each = nrow(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  # Maintain clustering order via factor levels
  plot_df$taxon <- factor(plot_df$taxon, levels = rev(rownames(mat)))
  plot_df$func <- factor(plot_df$func, levels = colnames(mat))

  text_size <- calculate_smart_text_size(max(nrow(mat), ncol(mat)),
                                         base_size = font_size)

  # Build heatmap
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$func, y = .data$taxon, fill = .data$value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient(
      low = low_color, high = high_color,
      name = "Mean\ncontribution"
    ) +
    ggplot2::labs(
      x = "Function",
      y = "Taxon",
      title = custom_title
    ) +
    ggplot2::theme_minimal(base_size = font_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, size = text_size
      ),
      axis.text.y = ggplot2::element_text(size = text_size),
      panel.grid = ggplot2::element_blank()
    )

  # Add dendrograms with patchwork if clustering is enabled
  if (!is.null(row_hc) || !is.null(col_hc)) {
    if (!is.null(row_hc)) {
      row_dendro <- create_dendrogram(row_hc, dendro_line_size = dendro_line_size,
                                      horizontal = TRUE)
    }
    if (!is.null(col_hc)) {
      col_dendro <- create_dendrogram(col_hc, dendro_line_size = dendro_line_size,
                                      horizontal = FALSE)
    }

    if (!is.null(row_hc) && !is.null(col_hc)) {
      # Both dendrograms
      if (!is.null(row_dendro) && !is.null(col_dendro)) {
        p <- (patchwork::plot_spacer() + col_dendro +
                row_dendro + p) +
          patchwork::plot_layout(
            ncol = 2,
            widths = c(0.15, 1),
            heights = c(0.15, 1)
          )
      } else {
        # Fallback if ggdendro not available
      }
    } else if (!is.null(row_hc) && !is.null(row_dendro)) {
      p <- (row_dendro + p) +
        patchwork::plot_layout(widths = c(0.15, 1))
    } else if (!is.null(col_hc) && !is.null(col_dendro)) {
      p <- (col_dendro / p) +
        patchwork::plot_layout(heights = c(0.15, 1))
    }
  }

  p
}
