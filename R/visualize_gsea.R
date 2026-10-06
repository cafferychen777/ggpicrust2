gsea_score_label <- function(gsea_results) {
  if ("score_label" %in% colnames(gsea_results)) {
    score_labels <- unique(as.character(gsea_results$score_label))
    score_labels <- score_labels[!is.na(score_labels) & nzchar(score_labels)]
    if (length(score_labels) == 1) {
      return(score_labels)
    }
    if (length(score_labels) > 1) {
      stop("gsea_results contains multiple score labels: ",
           paste(score_labels, collapse = ", "),
           ". Filter to one GSEA method or score type before plotting NES-based visualizations.",
           call. = FALSE)
    }
  }

  if ("score_type" %in% colnames(gsea_results)) {
    score_types <- unique(as.character(gsea_results$score_type))
    score_types <- score_types[!is.na(score_types) & nzchar(score_types)]
    if (length(score_types) > 1) {
      stop("gsea_results contains multiple score types: ",
           paste(score_types, collapse = ", "),
           ". Filter to one GSEA method or score type before plotting NES-based visualizations.",
           call. = FALSE)
    }
    if (length(score_types) == 1 && score_types == "signed_log10_pvalue") {
      return("Signed -log10(p-value)")
    }
  }

  "Normalized Enrichment Score (NES)"
}

#' Visualize GSEA results
#'
#' This function creates various visualizations for Gene Set Enrichment Analysis (GSEA) results.
#' It automatically detects whether pathway names are available (from gsea_pathway_annotation())
#' and uses them for better readability, falling back to pathway IDs if names are not available.
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param plot_type Visualization type: "enrichment_plot" (a score-summary
#' bar chart, not a running enrichment curve), "dotplot", "barplot",
#' "network", or "heatmap". Network edges use leading-edge overlap; the
#' heatmap shows per-pathway mean leading-edge abundance, standardized across
#' samples. These two modes need leading-edge genes from preranked methods,
#' which camera/fry do not provide.
#' @param n_pathways An integer specifying the number of pathways to display
#' @param sort_by A character string specifying the sorting criterion: "NES", "pvalue", or "p.adjust"
#' @param colors A vector of valid R colors used for the default heatmap group
#'   annotation palette. Colors are cycled when the number of observed groups
#'   exceeds the palette length.
#' @param abundance A data frame containing the original abundance data (required for heatmap visualization). Data frames may also provide a leading non-numeric feature ID column (for example \code{#NAME}, \code{feature}, or \code{pathway}); it is converted to row names before sample alignment.
#' @param metadata A data frame containing sample metadata (required for heatmap visualization)
#' @param group A character string specifying the column name in metadata that contains the grouping variable (required for heatmap visualization)
#' @param network_params A named list of network overrides. Supported names are
#'   `similarity_measure`, `similarity_cutoff`, `layout`, `node_color_by`, and
#'   `edge_width_by`.
#' @param heatmap_params A named list of heatmap overrides. Supported names are
#'   `cluster_rows`, `cluster_columns`, `show_rownames`, `annotation_colors`, and
#'   `col_fun`. Custom `annotation_colors` must be a list containing a named `Group`
#'   color vector whose names exactly match the observed group labels.
#' @param pathway_label_column A character string specifying which column to use for pathway labels.
#'   If NULL (default), the function will automatically use 'pathway_name' if available, otherwise 'pathway_id'.
#'   This allows for custom labeling when using annotated GSEA results.
#' @param scale Optional palette/scale for customizing colors. Accepts: (1) a character vector of colors,
#'   (2) a function that returns colors given an integer (e.g., viridisLite::viridis), or
#'   (3) for non-heatmap plots, a ggplot2 scale object for the mapped aesthetic
#'   (e.g., ggplot2::scale_fill_gradientn(...)).
#'   When NULL, defaults keep current behavior. Applies to: enrichment_plot (fill, continuous),
#'   dotplot (color, continuous), barplot (fill, discrete Positive/Negative), network (color, diverging around 0),
#'   heatmap (main heatmap col; row annotation stays default unless overridden in heatmap_params).
#'
#' For every \code{plot_type}, the selected rows after sorting and
#' \code{n_pathways} filtering must contain non-empty, unique
#' \code{pathway_id} values so that each displayed pathway maps to exactly
#' one GSEA result row.
#'
#' Results from \code{pathway_gsea(method = "camera")} and
#' \code{pathway_gsea(method = "fry")} include a legacy \code{NES} column for
#' visualization compatibility, but limma does not estimate a true normalized
#' enrichment score for these methods. When \code{score_label} is present,
#' axis and legend labels use it instead of labeling the value as NES.
#'
#' The function selects the top rows after sorting; it does not automatically
#' filter by significance. Inspect adjusted p-values before interpreting a
#' displayed pathway as significant.
#'
#' @return A ggplot2 object or ComplexHeatmap object
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#'
#' # Prepare abundance data
#' abundance_data <- as.data.frame(ko_abundance)
#' rownames(abundance_data) <- abundance_data[, "#NAME"]
#' abundance_data <- abundance_data[, -1]
#'
#' # Run GSEA analysis (using camera method - recommended)
#' gsea_results <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "camera"
#' )
#'
#' # Create enrichment plot with pathway IDs (default)
#' visualize_gsea(gsea_results, plot_type = "enrichment_plot", n_pathways = 10)
#'
#' # Annotate results for better pathway names
#' annotated_results <- gsea_pathway_annotation(
#'   gsea_results = gsea_results,
#'   pathway_type = "KEGG"
#' )
#'
#' # Create plots with readable pathway names
#' visualize_gsea(annotated_results, plot_type = "dotplot", n_pathways = 20)
#' visualize_gsea(annotated_results, plot_type = "barplot", n_pathways = 15)
#'
#' # Only preranked methods provide leading-edge genes for network/heatmap.
#' fgsea_results <- pathway_gsea(
#'   abundance_data, metadata, "Environment", method = "fgsea",
#'   comparison = c("Pro-survival", "Pro-inflammatory"), seed = 42
#' )
#' visualize_gsea(fgsea_results, plot_type = "network", n_pathways = 15)
#'
#' # Use custom column for labels (if available)
#' visualize_gsea(annotated_results, plot_type = "barplot",
#'                pathway_label_column = "pathway_name", n_pathways = 10)
#'
#' # Create heatmap
#' visualize_gsea(
#'   fgsea_results,
#'   plot_type = "heatmap",
#'   n_pathways = 15,
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment"
#' )
#' }
visualize_gsea <- function(gsea_results,
                          plot_type = "enrichment_plot",
                          n_pathways = 20,
                          sort_by = "p.adjust",
                          colors = NULL,
                          abundance = NULL,
                          metadata = NULL,
                          group = NULL,
                          network_params = list(),
                          heatmap_params = list(),
                          pathway_label_column = NULL,
                          scale = NULL) {

  # Input validation using unified functions
  validate_dataframe(gsea_results, param_name = "gsea_results")
  validate_choice(plot_type, c("enrichment_plot", "dotplot", "barplot", "network", "heatmap"), "plot_type")
  validate_choice(sort_by, c("NES", "pvalue", "p.adjust"), "sort_by")

  required_plot_cols <- switch(
    plot_type,
    enrichment_plot = c("pathway_id", "NES", "pvalue", "p.adjust"),
    dotplot = c("pathway_id", "NES", "p.adjust", "size"),
    barplot = c("pathway_id", "NES"),
    network = c("pathway_id", "NES", "pvalue", "p.adjust", "size", "leading_edge"),
    heatmap = c("pathway_id", "NES", "leading_edge")
  )
  required_cols <- unique(c(sort_by, required_plot_cols))
  validate_dataframe(gsea_results, required_cols = required_cols,
                     param_name = "gsea_results")
  validate_gsea_visualization_values(gsea_results, required_cols)
  score_label <- if ("NES" %in% required_cols) {
    gsea_score_label(gsea_results)
  } else {
    "Normalized Enrichment Score (NES)"
  }

  validate_color_values(colors, "colors", allow_null = TRUE)
  scale <- normalize_gsea_scale(scale, plot_type)
  if (!is.null(pathway_label_column) &&
      (!is.character(pathway_label_column) ||
       length(pathway_label_column) != 1 ||
       is.na(pathway_label_column) ||
       !nzchar(trimws(pathway_label_column)))) {
    stop("pathway_label_column must be NULL or a single non-empty character string",
         call. = FALSE)
  }
  validate_count_parameter(n_pathways, "n_pathways")
  n_pathways <- as.integer(n_pathways)

  gsea_results$pathway_id <- validate_nonempty_character_column(
    gsea_results$pathway_id,
    "pathway_id",
    "gsea_results"
  )

  # Note: enrichment_plot / dotplot / barplot are built entirely with
  # ggplot2 (see the branches below); they do not call into enrichplot
  # at runtime, so we don't require that Bioconductor package here.
  # Network and heatmap branches still depend on their own stacks and
  # are checked below.
  if (plot_type == "network") {
    network_params <- merge_named_parameters(
      defaults = list(
        similarity_measure = "jaccard",
        similarity_cutoff = 0.3,
        layout = "fruchterman",
        node_color_by = "NES",
        edge_width_by = "similarity"
      ),
      overrides = network_params,
      param_name = "network_params"
    )
    validate_choice(network_params$similarity_measure,
                    c("jaccard", "overlap", "correlation"),
                    "network_params$similarity_measure")
    validate_choice(network_params$layout,
                    c("fruchterman", "kamada", "circle"),
                    "network_params$layout")
    validate_choice(network_params$node_color_by,
                    c("NES", "pvalue", "p.adjust"),
                    "network_params$node_color_by")
    validate_choice(network_params$edge_width_by,
                    c("similarity", "constant"),
                    "network_params$edge_width_by")
    validate_probability_threshold(
      network_params$similarity_cutoff,
      "network_params$similarity_cutoff",
      allow_zero = TRUE
    )

    require_package("igraph", "network plots")
    require_package("ggraph", "network plots")
    require_package("tidygraph", "network plots")
  }

  if (plot_type == "heatmap") {
    if (is.null(abundance) || is.null(metadata) || is.null(group)) {
      stop("For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required",
           call. = FALSE)
    }

    heatmap_params <- merge_named_parameters(
      defaults = list(
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_rownames = TRUE,
        annotation_colors = NULL,
        col_fun = NULL
      ),
      overrides = heatmap_params,
      param_name = "heatmap_params"
    )
    heatmap_params$cluster_rows <- normalize_logical_flag(
      heatmap_params$cluster_rows,
      "heatmap_params$cluster_rows"
    )
    heatmap_params$cluster_columns <- normalize_logical_flag(
      heatmap_params$cluster_columns,
      "heatmap_params$cluster_columns"
    )
    heatmap_params$show_rownames <- normalize_logical_flag(
      heatmap_params$show_rownames,
      "heatmap_params$show_rownames"
    )
    if (!is.null(heatmap_params$col_fun) &&
        !is.function(heatmap_params$col_fun)) {
      stop("'heatmap_params$col_fun' must be NULL or a color function.",
           call. = FALSE)
    }
    if (!is.null(scale) && !is.null(heatmap_params$col_fun)) {
      stop(
        "Specify only one heatmap color mapping: 'scale' or 'heatmap_params$col_fun'.",
        call. = FALSE
      )
    }

    require_package("ComplexHeatmap", "heatmap plots")
    require_package("circlize", "heatmap plots")
  }

  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  }

  # Determine which column to use for pathway labels
  if (!is.null(pathway_label_column)) {
    # User specified a custom column
    require_column(gsea_results, pathway_label_column, "gsea_results")
    pathway_label_col <- pathway_label_column
  } else {
    # Auto-detect: prefer pathway_name if available, otherwise use pathway_id
    if ("pathway_name" %in% colnames(gsea_results)) {
      pathway_label_col <- "pathway_name"
    } else if ("pathway_id" %in% colnames(gsea_results)) {
      pathway_label_col <- "pathway_id"
    } else {
      # Check if we have any results first
      if (nrow(gsea_results) == 0) {
        return(create_empty_plot(plot_type))
      }
      stop("GSEA results must contain either 'pathway_name' or 'pathway_id' column")
    }
  }

  # Create a standardized pathway_label column for consistent use throughout
  # the function. Missing labels fall back row-wise to the stable pathway ID.
  pathway_labels <- as.character(gsea_results[[pathway_label_col]])
  missing_labels <- is.na(pathway_labels) | !nzchar(trimws(pathway_labels))
  pathway_labels[missing_labels] <- gsea_results$pathway_id[missing_labels]
  gsea_results$pathway_label <- pathway_labels

  # Sort results based on the specified criterion
  if (sort_by == "NES") {
    gsea_results <- gsea_results[order(abs(gsea_results$NES), decreasing = TRUE), ]
  } else if (sort_by == "pvalue") {
    gsea_results <- gsea_results[order(gsea_results$pvalue), ]
  } else if (sort_by == "p.adjust") {
    gsea_results <- gsea_results[order(gsea_results$p.adjust), ]
  }

  # Check if we have any results after filtering
  if (nrow(gsea_results) == 0) {
    return(create_empty_plot(plot_type))
  }

  # Limit to top n_pathways
  if (nrow(gsea_results) > n_pathways) {
    gsea_results <- head(gsea_results, n_pathways)
  }

  duplicated_pathways <- unique(gsea_results$pathway_id[duplicated(gsea_results$pathway_id)])
  if (length(duplicated_pathways) > 0) {
    stop(
      "visualize_gsea() requires one row per pathway_id after sorting and n_pathways filtering. ",
      "Duplicate pathway_id values: ",
      paste(utils::head(duplicated_pathways, 5), collapse = ", "),
      ". Filter to one GSEA method/contrast or deduplicate pathways before plotting.",
      call. = FALSE
    )
  }

  duplicated_labels <- duplicated(gsea_results$pathway_label) |
    duplicated(gsea_results$pathway_label, fromLast = TRUE)
  gsea_results$pathway_label[duplicated_labels] <- paste0(
    gsea_results$pathway_label[duplicated_labels],
    " [",
    gsea_results$pathway_id[duplicated_labels],
    "]"
  )
  if (anyDuplicated(gsea_results$pathway_label)) {
    # A literal label can itself contain an ID-like suffix. Falling back to the
    # unique primary key is preferable to inventing unstable positional labels.
    gsea_results$pathway_label <- gsea_results$pathway_id
  }

  # Create visualization based on plot_type
  if (plot_type == "enrichment_plot") {
    # Create a basic barplot of NES values.
    # Visual ordering is handled by reorder() in the aesthetic, so there is
    # no need to sort the data frame here. The user's sort_by parameter
    # already determined which pathways were selected (top N); the display
    # order within the plot is always by NES magnitude via the reorder() call.
    p <- ggplot2::ggplot(gsea_results, ggplot2::aes(x = reorder(.data$pathway_label, .data$NES), y = .data$NES, fill = .data$p.adjust)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      # apply user-provided scale if available (fallback to default)
      {
        sc <- .build_continuous_scale(aes = "fill", scale = scale, diverging = FALSE, name = "Adjusted p-value")
        if (is.null(sc)) ggplot2::scale_fill_gradient(low = "red", high = "blue") else sc
      } +
      ggplot2::labs(
        title = "GSEA Enrichment Results",
        x = "Pathway",
        y = score_label,
        fill = "Adjusted p-value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

  } else if (plot_type == "dotplot") {
    # Create dotplot
    # Sort by NES
    gsea_results$pathway_label <- factor(gsea_results$pathway_label,
                                      levels = gsea_results$pathway_label[order(gsea_results$NES)])

    p <- ggplot2::ggplot(gsea_results,
                       ggplot2::aes(x = .data$NES, y = .data$pathway_label, color = .data$p.adjust, size = .data$size)) +
      ggplot2::geom_point() +
      # apply user-provided scale if available (fallback to default)
      {
        sc <- .build_continuous_scale(aes = "color", scale = scale, diverging = FALSE, name = "Adjusted p-value")
        if (is.null(sc)) ggplot2::scale_color_gradient(low = "red", high = "blue") else sc
      } +
      ggplot2::labs(
        title = "GSEA Results",
        x = score_label,
        y = "Pathway",
        color = "Adjusted p-value",
        size = "Gene Set Size"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

  } else if (plot_type == "barplot") {
    # Create barplot
    # Sort by NES
    gsea_results$pathway_label <- factor(gsea_results$pathway_label,
                                      levels = gsea_results$pathway_label[order(gsea_results$NES)])

    # Add color based on NES direction
    gsea_results$direction <- ifelse(gsea_results$NES > 0, "Positive", "Negative")

    p <- ggplot2::ggplot(gsea_results,
                       ggplot2::aes(x = .data$pathway_label, y = .data$NES, fill = .data$direction)) +
      ggplot2::geom_bar(stat = "identity") +
      {
        sc <- .build_discrete_fill_for_direction(scale = scale)
        if (is.null(sc)) ggplot2::scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) else sc
      } +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "GSEA Results",
        x = "Pathway",
        y = score_label,
        fill = "Direction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

  } else if (plot_type == "network") {
    # Create network plot
    p <- create_network_plot(
      gsea_results = gsea_results,
      similarity_measure = network_params$similarity_measure,
      similarity_cutoff = network_params$similarity_cutoff,
      layout = network_params$layout,
      node_color_by = network_params$node_color_by,
      edge_width_by = network_params$edge_width_by,
      scale = scale
    )

  } else if (plot_type == "heatmap") {
    # Create heatmap
    p <- create_heatmap_plot(
      gsea_results = gsea_results,
      abundance = abundance,
      metadata = metadata,
      group = group,
      cluster_rows = heatmap_params$cluster_rows,
      cluster_columns = heatmap_params$cluster_columns,
      show_rownames = heatmap_params$show_rownames,
      annotation_colors = heatmap_params$annotation_colors,
      default_group_colors = colors,
      col_fun = {
        # Prefer explicit col_fun if provided in heatmap_params; else build from `scale`
        if (!is.null(heatmap_params$col_fun)) heatmap_params$col_fun else .build_heatmap_col_fun(scale)
      }
    )
  }

  return(p)
}

validate_gsea_visualization_values <- function(gsea_results, required_cols) {
  if ("NES" %in% required_cols) {
    validate_finite_numeric_values(gsea_results$NES, "NES", "gsea_results",
                                   allow_na = FALSE)
  }
  if ("pvalue" %in% required_cols) {
    validate_probability_values(gsea_results$pvalue, "pvalue", "gsea_results",
                                allow_na = FALSE)
  }
  if ("p.adjust" %in% required_cols) {
    validate_probability_values(gsea_results$p.adjust, "p.adjust", "gsea_results",
                                allow_na = FALSE)
  }
  if ("size" %in% required_cols) {
    validate_finite_numeric_values(gsea_results$size, "size", "gsea_results",
                                   allow_na = FALSE)
    bad_size <- gsea_results$size <= 0 |
      gsea_results$size > .Machine$integer.max |
      gsea_results$size != floor(gsea_results$size)
    if (any(bad_size)) {
      stop("Column 'size' in gsea_results must contain positive integer values.",
           call. = FALSE)
    }
  }
  if ("leading_edge" %in% required_cols) {
    if (is.factor(gsea_results$leading_edge)) {
      gsea_results$leading_edge <- as.character(gsea_results$leading_edge)
    }
    if (!is.character(gsea_results$leading_edge) &&
        !all(is.na(gsea_results$leading_edge))) {
      stop("Column 'leading_edge' in gsea_results must be a character vector or NA.",
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

parse_leading_edges <- function(leading_edge) {
  if (is.factor(leading_edge)) {
    leading_edge <- as.character(leading_edge)
  }
  leading_edge <- as.character(leading_edge)

  lapply(leading_edge, function(x) {
    if (length(x) != 1 || is.na(x) || !nzchar(trimws(x))) {
      return(character(0))
    }
    genes <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
    unique(genes[nzchar(genes) & !is.na(genes)])
  })
}

#' Create empty plot for edge cases
#'
#' @param plot_type A character string specifying the visualization type
#'
#' @return A ggplot2 object
#' @keywords internal
create_empty_plot <- function(plot_type) {
  # Create consistent empty plots for each visualization type
  base_plot <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = paste("No pathways to display for", plot_type))

  if (plot_type %in% c("network", "heatmap")) {
    # Special handling for complex plot types
    base_plot + ggplot2::annotate("text", x = 0, y = 0,
                                  label = "No significant pathways found",
                                  size = 4, color = "gray50")
  } else {
    base_plot + ggplot2::annotate("text", x = 0, y = 0,
                                  label = "No pathways to display",
                                  size = 4, color = "gray50")
  }
}

#' Internal: detect if an object is a ggplot2 Scale
#' @keywords internal
.is_ggplot_scale <- function(x) {
  inherits(x, "Scale")
}

#' Internal: validate and normalize the GSEA visualization scale
#' @noRd
normalize_gsea_scale <- function(scale, plot_type) {
  if (is.null(scale)) {
    return(NULL)
  }

  if (.is_ggplot_scale(scale)) {
    if (identical(plot_type, "heatmap")) {
      stop(
        "A ggplot2 scale object cannot be used for a ComplexHeatmap plot. Supply a color vector, palette function, or 'heatmap_params$col_fun'.",
        call. = FALSE
      )
    }

    required_aesthetic <- if (plot_type %in% c("enrichment_plot", "barplot")) {
      "fill"
    } else {
      "colour"
    }
    scale_aesthetics <- as.character(scale$aesthetics)
    if (required_aesthetic == "colour") {
      scale_aesthetics[scale_aesthetics == "color"] <- "colour"
    }
    if (!required_aesthetic %in% scale_aesthetics) {
      stop(
        "The ggplot2 object supplied as 'scale' must map the '",
        required_aesthetic,
        "' aesthetic for plot_type = '",
        plot_type,
        "'.",
        call. = FALSE
      )
    }
    if (identical(plot_type, "barplot") &&
        !inherits(scale, "ScaleDiscrete")) {
      stop("The ggplot2 scale for a barplot must be discrete.", call. = FALSE)
    }
    if (!identical(plot_type, "barplot") &&
        inherits(scale, "ScaleDiscrete")) {
      stop(
        "The ggplot2 scale for '", plot_type,
        "' must be continuous because it maps numeric values.",
        call. = FALSE
      )
    }
    return(scale)
  }

  if (is.function(scale)) {
    palette_result <- tryCatch(
      list(colors = scale(100L), error = NULL),
      error = function(e) list(colors = NULL, error = e)
    )
    if (!is.null(palette_result$error)) {
      stop(
        "The palette function supplied as 'scale' failed: ",
        conditionMessage(palette_result$error),
        call. = FALSE
      )
    }
    scale <- palette_result$colors
  }

  if (!is.character(scale)) {
    stop(
      "'scale' must be NULL, a character color vector, a palette function, or a compatible ggplot2 scale object.",
      call. = FALSE
    )
  }
  validate_color_values(scale, "scale")
  if (length(scale) < 2) {
    stop("'scale' must provide at least two colors.", call. = FALSE)
  }

  scale
}

#' Internal: extract a normalized color vector
#' @keywords internal
.as_color_vector <- function(scale) {
  if (is.null(scale)) {
    return(NULL)
  }
  if (!is.character(scale)) {
    stop("Internal error: color scale was not normalized.", call. = FALSE)
  }
  scale
}

#' Internal: build a continuous ggplot2 scale layer from colors or ggplot2 scale
#' @keywords internal
.build_continuous_scale <- function(aes = c("fill", "color"), scale = NULL,
                                   diverging = FALSE, midpoint = NULL, name = NULL) {
  aes <- match.arg(aes)
  # If a ggplot2 scale object is provided, return as is
  if (.is_ggplot_scale(scale)) return(scale)
  cols <- .as_color_vector(scale)
  if (is.null(cols)) return(NULL)
  if (diverging && length(cols) >= 3) {
    low <- cols[1]
    mid <- cols[ceiling(length(cols) / 2)]
    high <- cols[length(cols)]
    if (aes == "fill") {
      return(ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high,
                                           midpoint = if (is.null(midpoint)) 0 else midpoint,
                                           name = name))
    } else {
      return(ggplot2::scale_color_gradient2(low = low, mid = mid, high = high,
                                            midpoint = if (is.null(midpoint)) 0 else midpoint,
                                            name = name))
    }
  } else {
    # General gradientn for sequential palettes
    if (aes == "fill") {
      return(ggplot2::scale_fill_gradientn(colors = cols, name = name))
    } else {
      return(ggplot2::scale_color_gradientn(colors = cols, name = name))
    }
  }
}

#' Internal: build a discrete fill scale for barplot direction
#' @keywords internal
.build_discrete_fill_for_direction <- function(scale = NULL) {
  if (.is_ggplot_scale(scale)) return(scale)
  cols <- .as_color_vector(scale)
  if (is.null(cols)) return(NULL)
  values <- c("Positive" = cols[length(cols)], "Negative" = cols[1])
  ggplot2::scale_fill_manual(values = values)
}

#' Internal: build a circlize colorRamp2 function for ComplexHeatmap from user scale
#' @keywords internal
.build_heatmap_col_fun <- function(scale = NULL) {
  # The public boundary has already normalized functions to color vectors and
  # rejected ggplot2 scales for ComplexHeatmap output.
  cols <- .as_color_vector(scale)

  # If no colors provided or couldn't convert, return NULL (use defaults)
  if (is.null(cols)) return(NULL)

  # Create a diverging color function for heatmap
  # Typical heatmap shows z-scores, so we use -2, 0, 2 as breakpoints
  if (length(cols) >= 3) {
    # Use low, mid, high colors from the palette
    low_col <- cols[1]
    mid_col <- cols[ceiling(length(cols) / 2)]
    high_col <- cols[length(cols)]
    return(circlize::colorRamp2(c(-2, 0, 2), c(low_col, mid_col, high_col)))
  } else if (length(cols) == 2) {
    # If only 2 colors, create a gradient without midpoint
    return(circlize::colorRamp2(c(-2, 2), c(cols[1], cols[2])))
  }
}

#' Create network visualization of GSEA results
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param similarity_measure A character string specifying the similarity measure: "jaccard", "overlap", or "correlation"
#' @param similarity_cutoff A numeric value specifying the similarity threshold for filtering connections
#' @param layout A character string specifying the network layout algorithm: "fruchterman", "kamada", or "circle"
#' @param node_color_by A character string specifying the node color mapping: "NES", "pvalue", or "p.adjust"
#' @param edge_width_by A character string specifying the edge width mapping: "similarity" or "constant"
#' @param scale Optional palette/scale for customizing node color mapping (same conventions as visualize_gsea)
#'
#' @return A ggplot2 object
#' @keywords internal
create_network_plot <- function(gsea_results,
                               similarity_measure = "jaccard",
                               similarity_cutoff = 0.3,
                               layout = "fruchterman",
                               node_color_by = "NES",
                               edge_width_by = "similarity",
                               scale = NULL) {
  # Note: Input validation (packages, n_pathways, nrow) is done in visualize_gsea()
  node_color_label <- if (identical(node_color_by, "NES")) {
    gsea_score_label(gsea_results)
  } else {
    node_color_by
  }

  # Extract leading edge genes. Missing values are unknown, not a literal shared
  # gene token; parse them as empty sets so they cannot create false network
  # similarities between pathways with missing leading-edge information.
  leading_edges <- parse_leading_edges(gsea_results$leading_edge)
  names(leading_edges) <- gsea_results$pathway_id

  # Calculate pathway similarity
  n <- length(leading_edges)
  pathway_ids <- names(leading_edges)
  similarity_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(similarity_matrix) <- pathway_ids
  colnames(similarity_matrix) <- pathway_ids

  for (i in seq_len(max(n - 1, 0))) {
    for (j in seq.int(i + 1, n)) {
      set1 <- leading_edges[[i]]
      set2 <- leading_edges[[j]]

      if (length(set1) == 0 || length(set2) == 0) {
        next
      }

      intersection_size <- length(intersect(set1, set2))
      similarity <- if (similarity_measure == "jaccard") {
        intersection_size / length(union(set1, set2))
      } else if (similarity_measure == "overlap") {
        intersection_size / min(length(set1), length(set2))
      } else {
        # Binary cosine similarity, retained under the legacy "correlation"
        # option for backward compatibility.
        intersection_size / sqrt(length(set1) * length(set2))
      }
      similarity_matrix[i, j] <- similarity
      similarity_matrix[j, i] <- similarity
    }
  }

  # Apply similarity cutoff
  similarity_matrix[similarity_matrix < similarity_cutoff] <- 0

  # Check if there are any connections after applying cutoff
  if (!any(similarity_matrix > 0)) {
    return(create_empty_plot("network"))
  }

  # Create graph object
  graph <- igraph::graph_from_adjacency_matrix(
    similarity_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # Add node attributes
  vertex_attr <- data.frame(
    name = pathway_ids,
    NES = gsea_results$NES[match(pathway_ids, gsea_results$pathway_id)],
    pvalue = gsea_results$pvalue[match(pathway_ids, gsea_results$pathway_id)],
    p.adjust = gsea_results$p.adjust[match(pathway_ids, gsea_results$pathway_id)],
    size = gsea_results$size[match(pathway_ids, gsea_results$pathway_id)],
    pathway_label = gsea_results$pathway_label[match(pathway_ids, gsea_results$pathway_id)],
    stringsAsFactors = FALSE
  )

  # Create a tidygraph object
  tbl_graph <- tidygraph::as_tbl_graph(graph) %>%
    tidygraph::activate("nodes") %>%
    dplyr::mutate(
      name = vertex_attr$name,
      NES = vertex_attr$NES,
      pvalue = vertex_attr$pvalue,
      p.adjust = vertex_attr$p.adjust,
      size = vertex_attr$size,
      pathway_label = vertex_attr$pathway_label
    )

  # Select layout algorithm
  if (layout == "fruchterman") {
    layout_name <- "fr"
  } else if (layout == "kamada") {
    layout_name <- "kk"
  } else if (layout == "circle") {
    layout_name <- "circle"
  } else {
    layout_name <- "fr"
  }

  edge_layer <- if (edge_width_by == "similarity") {
    ggraph::geom_edge_link(ggplot2::aes(width = .data$weight, alpha = .data$weight))
  } else {
    ggraph::geom_edge_link(ggplot2::aes(alpha = .data$weight), width = 0.5)
  }
  edge_width_scale <- if (edge_width_by == "similarity") {
    ggraph::scale_edge_width(range = c(0.1, 2))
  } else {
    NULL
  }

  # Create ggraph visualization
  p <- ggraph::ggraph(tbl_graph, layout = layout_name) +
    edge_layer +
    ggraph::geom_node_point(ggplot2::aes(color = .data[[node_color_by]], size = .data$size)) +
    ggraph::geom_node_text(ggplot2::aes(label = .data$pathway_label), repel = TRUE, size = 3) +
    edge_width_scale +
    ggraph::scale_edge_alpha(range = c(0.1, 0.8)) +
    # apply user-provided diverging scale for node color if available (fallback to default)
    {
      sc <- .build_continuous_scale(aes = "color", scale = scale, diverging = TRUE, midpoint = 0, name = node_color_label)
      if (is.null(sc)) ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = node_color_label) else sc
    } +
    ggplot2::scale_size(range = c(2, 8), name = "Gene Set Size") +
    ggraph::theme_graph() +
    ggplot2::labs(
      title = "GSEA Pathway Network",
      subtitle = paste("Similarity measure:", similarity_measure, "| Cutoff:", similarity_cutoff)
    )

  return(p)
}

#' Create heatmap visualization of GSEA results
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param abundance A data frame containing the original abundance data
#' @param metadata A data frame containing sample metadata
#' @param group A character string specifying the column name in metadata that contains the grouping variable
#' @param cluster_rows A logical value indicating whether to cluster rows
#' @param cluster_columns A logical value indicating whether to cluster columns
#' @param show_rownames A logical value indicating whether to show row names
#' @param annotation_colors A list of colors for annotations
#' @param default_group_colors Colors used to construct the group annotation
#'   palette when annotation_colors is NULL
#' @param col_fun A color function (e.g., circlize::colorRamp2) to control the main heatmap colors (optional)
#'
#' @return A ComplexHeatmap object
#' @keywords internal
create_heatmap_plot <- function(gsea_results,
                               abundance,
                               metadata,
                               group,
                               cluster_rows = TRUE,
                               cluster_columns = TRUE,
                               show_rownames = TRUE,
                               annotation_colors = NULL,
                               default_group_colors = c("#E41A1C", "#377EB8"),
                               col_fun = NULL) {
  # Note: Input validation (packages, n_pathways, nrow) is done in visualize_gsea()
  score_label <- gsea_score_label(gsea_results)

  abundance <- normalize_abundance_feature_ids(
    abundance,
    context = "visualize_gsea() heatmap abundance"
  )

  # Extract leading edge genes
  leading_edges <- parse_leading_edges(gsea_results$leading_edge)
  names(leading_edges) <- gsea_results$pathway_id

  # Check if any leading edges are empty
  if (all(lengths(leading_edges) == 0)) {
    return(ComplexHeatmap::Heatmap(
      matrix(0, nrow = 1, ncol = 1),
      name = "Empty",
      show_row_names = FALSE,
      show_column_names = FALSE,
      row_title = "No leading edge genes found",
      column_title = "No gene expression data"
    ))
  }

  # Align abundance columns with metadata rows using the package-wide
  # sample-alignment utility. Previously this branch relied solely on
  # `rownames(metadata)` to locate samples, which silently broke on the
  # common case of metadata carrying a `sample_name`/`sample_id` column
  # with default integer rownames -- the column annotation then came
  # out as all-NA, diverging from how every other function in the
  # package matches abundance to metadata.
  aligned <- align_samples(abundance, metadata, verbose = FALSE)
  abundance <- as.matrix(aligned$abundance)
  validate_numeric_matrix(abundance, check_negative = FALSE)
  validate_finite_numeric_values(as.vector(abundance), "abundance", "abundance",
                                 allow_na = FALSE)
  metadata <- aligned$metadata
  validate_group(metadata, group, min_groups = 1)
  validate_group_vector_for_summary(
    metadata[[group]],
    context = paste0("Group column '", group, "' after sample alignment"),
    sample_ids = colnames(abundance),
    min_groups = 1
  )

  group_levels <- unique(as.character(metadata[[group]]))
  if (is.null(annotation_colors)) {
    group_palette <- rep(default_group_colors, length.out = length(group_levels))
    annotation_colors <- list(
      Group = normalize_level_colors(
        group_palette,
        group_levels,
        "heatmap group colors"
      )
    )
  } else {
    if (!is.list(annotation_colors) ||
        !identical(names(annotation_colors), "Group")) {
      stop(
        "'heatmap_params$annotation_colors' must be NULL or a list with exactly one named element, 'Group'.",
        call. = FALSE
      )
    }
    annotation_colors$Group <- normalize_level_colors(
      annotation_colors$Group,
      group_levels,
      "heatmap_params$annotation_colors$Group",
      allow_unnamed = FALSE
    )
  }

  # Create heatmap data matrix
  # For each pathway, calculate the average expression of leading edge genes
  heatmap_data <- matrix(0, nrow = length(leading_edges), ncol = ncol(abundance))
  rownames(heatmap_data) <- names(leading_edges)
  colnames(heatmap_data) <- colnames(abundance)

  matched_gene_counts <- integer(length(leading_edges))
  names(matched_gene_counts) <- names(leading_edges)

  for (i in seq_along(leading_edges)) {
    genes <- leading_edges[[i]]

    # Ensure all genes are in abundance data
    genes <- genes[genes %in% rownames(abundance)]
    matched_gene_counts[i] <- length(genes)

    if (length(genes) > 0) {
      # Calculate average abundance
      heatmap_data[i, ] <- colMeans(abundance[genes, , drop = FALSE])
    }
  }

  nonempty_sets <- lengths(leading_edges) > 0
  if (any(nonempty_sets) && all(matched_gene_counts[nonempty_sets] == 0)) {
    stop(
      "None of the non-empty leading-edge genes were found in abundance row names. ",
      "Check that abundance row names use the same gene/KO identifiers as gsea_results$leading_edge.",
      call. = FALSE
    )
  }

  missing_sets <- names(matched_gene_counts)[nonempty_sets & matched_gene_counts == 0]
  if (length(missing_sets) > 0) {
    warning(
      length(missing_sets), " pathway(s) had non-empty leading edges but no matching abundance rows: ",
      paste(head(missing_sets, 5), collapse = ", "),
      call. = FALSE
    )
  }

  # Scale data. Constant rows produce NA after t(scale(t(.))); coerce
  # them to 0 so ComplexHeatmap doesn't crash on clustering.
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  heatmap_data_scaled[!is.finite(heatmap_data_scaled)] <- 0

  # Column annotation now derives from the aligned metadata, which is
  # guaranteed to match `colnames(heatmap_data)` by construction.
  column_annotation <- metadata[[group]]
  names(column_annotation) <- colnames(heatmap_data)

  # Create column annotation object
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Group = column_annotation,
    col = list(Group = annotation_colors$Group),
    show_legend = TRUE
  )

  # Create row annotation (pathway enrichment scores)
  row_annotation <- data.frame(
    NES = gsea_results$NES,
    row.names = gsea_results$pathway_id
  )

  # Ensure row annotation matches heatmap rows
  row_annotation <- row_annotation[rownames(heatmap_data), , drop = FALSE]

  # Create row annotation object. Use a generic annotation name when the score
  # is not a true normalized enrichment score, but keep the legend title
  # explicit through the annotation label.
  # Use symmetric breaks centered at 0 to ensure colorRamp2 gets
  # strictly increasing values even when all NES are same sign
  nes_abs_max <- max(abs(row_annotation$NES), 0.1)
  score_annotation_name <- if (identical(score_label, "Normalized Enrichment Score (NES)")) {
    "NES"
  } else {
    "Score"
  }
  score_values <- stats::setNames(list(row_annotation$NES), score_annotation_name)
  score_colors <- stats::setNames(
    list(circlize::colorRamp2(
      c(-nes_abs_max, 0, nes_abs_max),
      c("blue", "white", "red")
    )),
    score_annotation_name
  )
  legend_params <- stats::setNames(
    list(list(title = score_label)),
    score_annotation_name
  )
  ra <- do.call(
    ComplexHeatmap::rowAnnotation,
    c(score_values, list(
    col = score_colors,
    annotation_legend_param = legend_params,
    show_legend = TRUE
    ))
  )

  # Create heatmap
  heatmap <- ComplexHeatmap::Heatmap(
    heatmap_data_scaled,
    name = "Z-score",
    col = {
      # Use user-provided col_fun if present; else default blue-white-red
      if (!is.null(col_fun)) col_fun else circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    },
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_rownames,
    row_names_gp = grid::gpar(fontsize = 8),
    top_annotation = ha,
    right_annotation = ra,
    row_title = "Pathways",
    column_title = "Samples",
    row_names_max_width = grid::unit(15, "cm")
  )

  return(heatmap)
}
