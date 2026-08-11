# =============================================================================
# Taxa-Function Contribution Visualization
# =============================================================================
# Stacked bar plots and heatmaps showing which taxa drive pathway abundances.

#' Complete a contribution grid with structural zeros
#'
#' @noRd
complete_contribution_grid <- function(complete_index, observed, key_cols) {
  if (anyDuplicated(observed[key_cols])) {
    stop("Internal error: observed contribution keys must be unique.",
         call. = FALSE)
  }
  completed <- dplyr::left_join(complete_index, observed, by = key_cols)
  completed$contribution[is.na(completed$contribution)] <- 0
  completed
}

#' Complete sample/function contribution totals with structural zeros
#'
#' @noRd
complete_sample_function_totals <- function(contrib_data, samples,
                                            function_ids) {
  observed <- stats::aggregate(
    contribution ~ sample + function_id,
    data = contrib_data,
    FUN = sum
  )
  complete_index <- expand.grid(
    sample = samples,
    function_id = function_ids,
    stringsAsFactors = FALSE
  )
  complete_contribution_grid(
    complete_index,
    observed,
    key_cols = c("sample", "function_id")
  )
}

#' Stacked bar plot of taxa contributions
#'
#' Creates a stacked bar plot showing taxa contributions to predicted
#' functional abundances, faceted by function or sample group.
#'
#' @param contrib_agg A data.frame from \code{\link{aggregate_taxa_contributions}}.
#' @param metadata A data.frame containing sample metadata.
#' @param group Character. Column name in \code{metadata} for grouping samples.
#' @param function_ids Optional character vector of function IDs to plot.
#'   If NULL (default), the top \code{n_functions} by between-sample variance
#'   in total contribution are shown. Single-sample inputs are ranked by total
#'   contribution because variance is undefined. Facets preserve this ranking
#'   or the order of explicitly supplied IDs.
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
#' @details
#' The \code{sample}, \code{function_id}, and \code{taxon_label} columns must
#' contain non-empty values without \code{NA}. These columns define plotting and
#' aggregation groups, so missing identifiers would otherwise be dropped by R
#' aggregation or shown as unlabeled categories.
#' When \code{show_percentage = TRUE}, every plotted sample/function
#' combination must have a positive total contribution. Relative percentages
#' are undefined for zero-total combinations; use \code{show_percentage = FALSE}
#' to display absolute zero contributions.
#'
#' @examples
#' \donttest{
#' # Synthetic example
#' agg <- expand.grid(
#'   sample = c("S1", "S2", "S3", "S4"),
#'   function_id = c("K00001", "K00002"),
#'   taxon_label = c("Genus_A", "Genus_B", "Other"),
#'   stringsAsFactors = FALSE
#' )
#' agg$contribution <- runif(nrow(agg))
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
  if (nrow(contrib_agg) == 0) {
    stop("'contrib_agg' must contain at least one contribution row.",
         call. = FALSE)
  }
  validate_contrib_key_columns(contrib_agg,
                               c("sample", "function_id", "taxon_label"),
                               "contrib_agg")
  validate_contribution_values(contrib_agg$contribution, "contribution")
  validate_metadata(metadata)
  validate_group(metadata, group, min_groups = 1)
  validate_count_parameter(n_functions, "n_functions")
  show_percentage <- normalize_logical_flag(show_percentage,
                                            "show_percentage")
  validate_positive_number(font_size, "font_size")
  validate_choice(
    legend_position,
    c("top", "bottom", "left", "right", "none"),
    "legend_position"
  )
  function_ids_requested <- !is.null(function_ids)
  if (!is.null(function_ids)) {
    function_ids <- unique(validate_nonempty_character_column(
      function_ids,
      "function_ids",
      "function_ids"
    ))
  }

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
  group_values <- validate_group_vector_for_summary(
    metadata_aligned[[group]],
    context = paste0("metadata column '", group, "' after sample alignment"),
    sample_ids = common_samples,
    min_groups = 1
  )

  contrib_agg <- contrib_agg[contrib_agg$sample %in% common_samples, ]
  if (nrow(contrib_agg) == 0) {
    stop("No overlapping samples between contrib_agg and metadata.")
  }

  # Select functions to plot
  if (is.null(function_ids)) {
    sample_function_totals <- complete_sample_function_totals(
      contrib_agg,
      samples = unique(contrib_agg$sample),
      function_ids = unique(contrib_agg$function_id)
    )

    sample_count <- length(unique(sample_function_totals$sample))
    ranking_fun <- if (sample_count > 1) stats::var else sum
    func_var <- stats::aggregate(
      contribution ~ function_id,
      data = sample_function_totals,
      FUN = ranking_fun
    )
    func_var <- func_var[order(-func_var$contribution), ]
    function_ids <- utils::head(func_var$function_id, n_functions)
  }
  if (function_ids_requested) {
    missing_function_ids <- setdiff(function_ids, unique(contrib_agg$function_id))
    if (length(missing_function_ids) > 0) {
      stop(
        "No contribution rows match the requested function_ids after sample ",
        "alignment: ",
        paste(utils::head(missing_function_ids, 5), collapse = ", "),
        call. = FALSE
      )
    }
  }
  contrib_agg <- contrib_agg[contrib_agg$function_id %in% function_ids, ]
  if (nrow(contrib_agg) == 0) {
    stop("No contribution rows match the requested function_ids.",
         call. = FALSE)
  }
  contrib_agg$function_id <- factor(
    contrib_agg$function_id,
    levels = function_ids
  )

  # Add group info
  sample_col <- aligned$sample_col
  group_map <- stats::setNames(
    group_values,
    metadata_aligned[[sample_col]]
  )
  contrib_agg$group_var <- group_map[contrib_agg$sample]

  # Normalize to percentage if requested
  if (show_percentage) {
    sample_function_totals <- complete_sample_function_totals(
      contrib_agg,
      samples = common_samples,
      function_ids = function_ids
    )
    zero_totals <- sample_function_totals$contribution <= 0
    if (any(zero_totals)) {
      zero_examples <- paste(
        sample_function_totals$sample[zero_totals],
        sample_function_totals$function_id[zero_totals],
        sep = "/"
      )
      stop(
        "Cannot compute relative contribution percentages for sample/function ",
        "combinations with total contribution <= 0: ",
        paste(utils::head(zero_examples, 5), collapse = ", "),
        ". Use show_percentage = FALSE or remove zero-total combinations.",
        call. = FALSE
      )
    }

    total_col <- ".ggpicrust2_contribution_total"
    while (total_col %in% colnames(contrib_agg)) {
      total_col <- paste0(total_col, "_")
    }
    normalization_totals <- sample_function_totals
    colnames(normalization_totals)[
      colnames(normalization_totals) == "contribution"
    ] <- total_col
    contrib_agg <- dplyr::left_join(
      contrib_agg,
      normalization_totals,
      by = c("sample", "function_id")
    )
    contrib_agg$contribution <-
      contrib_agg$contribution / contrib_agg[[total_col]] * 100
    contrib_agg[[total_col]] <- NULL
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
#'   for replacing function IDs with readable descriptions. It must contain
#'   either \code{feature}/\code{description} or
#'   \code{pathway}/\code{pathway_name} columns.
#' @param n_functions Integer. Number of functions to include. Default 20.
#' @param cluster_rows Logical. Cluster rows (taxa)? Default TRUE.
#' @param cluster_cols Logical. Cluster columns (functions)? Default TRUE.
#' @param clustering_method Character. Method for \code{hclust}. Default \code{"complete"}. Ward methods require \code{clustering_distance = "euclidean"}.
#' @param clustering_distance Character. Distance metric. Default \code{"euclidean"}. Supported values are \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"}, and \code{"minkowski"}.
#' @param low_color Character. Color for low values. Default \code{"#f7f7f7"}.
#' @param high_color Character. Color for high values. Default \code{"#ca0020"}.
#' @param font_size Numeric. Base font size. Default 12.
#' @param dendro_line_size Numeric. Dendrogram line width. Default 0.5.
#' @param custom_title Optional plot title.
#'
#' @return A \code{ggplot2} or \code{patchwork} object.
#'
#' @details
#' The \code{sample}, \code{function_id}, and \code{taxon_label} columns must
#' contain non-empty values without \code{NA}. These columns define plotting and
#' aggregation groups, so missing identifiers would otherwise be dropped by R
#' aggregation or shown as unlabeled categories.
#' PICRUSt2 contribution outputs are sparse; combinations absent from
#' \code{contrib_agg} are treated as zero when computing mean contribution
#' across samples.
#' If \code{annotation_data} contains multiple non-empty labels for the same
#' plotted function ID, the function errors instead of silently choosing one
#' label. Repeated rows with the same ID and same label are allowed, and label
#' whitespace is normalized before comparison and display.
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
  if (nrow(contrib_agg) == 0) {
    stop("'contrib_agg' must contain at least one contribution row.",
         call. = FALSE)
  }
  validate_contrib_key_columns(contrib_agg,
                               c("sample", "function_id", "taxon_label"),
                               "contrib_agg")
  validate_count_parameter(n_functions, "n_functions")
  validate_contribution_values(contrib_agg$contribution, "contribution")
  cluster_rows <- normalize_logical_flag(cluster_rows, "cluster_rows")
  cluster_cols <- normalize_logical_flag(cluster_cols, "cluster_cols")
  validate_positive_number(font_size, "font_size")
  validate_positive_number(dendro_line_size, "dendro_line_size",
                           allow_zero = TRUE)
  validate_color_values(low_color, "low_color", expected_length = 1)
  validate_color_values(high_color, "high_color", expected_length = 1)
  validate_hclust_parameters(
    clustering_method,
    clustering_distance,
    allow_correlation = FALSE
  )

  # Select top functions by total contribution
  func_totals <- stats::aggregate(
    contribution ~ function_id, data = contrib_agg, FUN = sum
  )
  func_totals <- func_totals[order(-func_totals$contribution), ]
  top_funcs <- utils::head(func_totals$function_id, n_functions)
  heatmap_data <- contrib_agg[contrib_agg$function_id %in% top_funcs, ]

  # PICRUSt2 contribution outputs are sparse: zero contribution combinations
  # are commonly absent rather than represented as explicit rows. Treat absent
  # sample/function/taxon combinations as zero before computing sample means.
  sample_contrib <- stats::aggregate(
    contribution ~ sample + function_id + taxon_label,
    data = heatmap_data,
    FUN = sum
  )
  complete_index <- expand.grid(
    sample = unique(contrib_agg$sample),
    function_id = top_funcs,
    taxon_label = unique(heatmap_data$taxon_label),
    stringsAsFactors = FALSE
  )
  complete_contrib <- complete_contribution_grid(
    complete_index,
    sample_contrib,
    key_cols = c("sample", "function_id", "taxon_label")
  )

  mean_contrib <- stats::aggregate(
    contribution ~ function_id + taxon_label,
    data = complete_contrib,
    FUN = mean
  )

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
  mat <- mat[, top_funcs[top_funcs %in% colnames(mat)], drop = FALSE]

  # Replace function IDs with annotations if available.
  # pathway_annotation() produces `feature` as the ID column and
  # `description` as the readable label (non-ko_to_kegg mode); with
  # ko_to_kegg=TRUE it produces `pathway_name`. Accept either shape so
  # users can drop in whatever pathway_annotation() returned without
  # silently getting raw IDs on the heatmap axis.
  if (!is.null(annotation_data)) {
    annotation_cols <- resolve_contribution_annotation_columns(annotation_data)
    desc_map <- build_annotation_label_map(
      annotation_data,
      id_col = annotation_cols$id,
      label_col = annotation_cols$label,
      selected_ids = colnames(mat)
    )
    function_ids <- colnames(mat)
    new_names <- desc_map[function_ids]
    # Only replace where we found a non-empty match, truncate long names
    found <- !is.na(new_names) & nzchar(new_names)
    new_names[found] <- substr(new_names[found], 1, 50)
    new_names[!found] <- function_ids[!found]

    # Different function IDs can legitimately share the same annotation.
    # Axis factor levels must still be unique, so retain the function ID for
    # every ambiguous display label rather than crashing in factor().
    duplicated_labels <- duplicated(new_names) |
      duplicated(new_names, fromLast = TRUE)
    new_names[duplicated_labels] <- paste0(
      substr(new_names[duplicated_labels], 1, 38),
      " [", function_ids[duplicated_labels], "]"
    )
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

#' Resolve a supported taxa-contribution annotation schema
#'
#' @noRd
resolve_contribution_annotation_columns <- function(annotation_data) {
  validate_dataframe(annotation_data, param_name = "annotation_data")
  schemas <- list(
    c(id = "feature", label = "description"),
    c(id = "pathway", label = "pathway_name")
  )
  matched <- vapply(
    schemas,
    function(schema) all(unname(schema) %in% colnames(annotation_data)),
    logical(1)
  )
  if (!any(matched)) {
    stop(
      "'annotation_data' must contain either 'feature'/'description' or ",
      "'pathway'/'pathway_name' columns.",
      call. = FALSE
    )
  }

  as.list(schemas[[which(matched)[1]]])
}

#' Build a unique annotation label map for selected function IDs
#'
#' @noRd
build_annotation_label_map <- function(annotation_data, id_col, label_col,
                                       selected_ids) {
  annotation_ids <- as.character(annotation_data[[id_col]])
  annotation_labels <- trimws(as.character(annotation_data[[label_col]]))
  selected <- !is.na(annotation_ids) & annotation_ids %in% selected_ids
  if (!any(selected)) {
    return(character(0))
  }

  mapping <- data.frame(
    id = annotation_ids[selected],
    label = annotation_labels[selected],
    stringsAsFactors = FALSE
  )
  valid_label <- !is.na(mapping$label) & nzchar(trimws(mapping$label))
  mapping <- unique(mapping[valid_label, , drop = FALSE])
  if (nrow(mapping) == 0) {
    return(character(0))
  }

  labels_by_id <- split(mapping$label, mapping$id)
  conflicting_ids <- names(labels_by_id)[vapply(
    labels_by_id,
    function(labels) length(unique(labels)) > 1,
    logical(1)
  )]
  if (length(conflicting_ids) > 0) {
    stop(
      "annotation_data contains multiple labels for the same function ID: ",
      paste(utils::head(conflicting_ids, 5), collapse = ", "),
      ". Provide one annotation label per function ID before plotting.",
      call. = FALSE
    )
  }

  vapply(labels_by_id, function(labels) unique(labels)[1], character(1))
}
