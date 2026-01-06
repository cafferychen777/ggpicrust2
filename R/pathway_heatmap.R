#' Create pathway heatmap with support for multiple grouping variables
#'
#' This function creates a heatmap of the predicted functional pathway abundance data
#' with support for single or multiple grouping variables. The function first performs
#' z-score normalization on the abundance data, then converts it to a long format and
#' orders the samples based on the grouping information. The heatmap supports nested
#' faceting for multiple grouping variables and is created using the `ggplot2` library.
#'
#' @name pathway_heatmap
#' @param abundance A matrix or data frame of pathway abundance data, where columns
#'   correspond to samples and rows correspond to pathways. Must contain at least
#'   two samples.
#' @param metadata A data frame of metadata, where each row corresponds to a sample
#'   and each column corresponds to a metadata variable.
#' @param group A character string specifying the column name in the metadata data frame
#'   that contains the primary group variable. Must contain at least two groups.
#' @param secondary_groups A character vector specifying additional grouping variables
#'   for creating nested faceted heatmaps. If NULL, only the primary group will be used.
#'   These variables will be used as secondary levels in the faceting hierarchy.
#' @param colors A vector of colors used for the background of the facet labels in the
#'   heatmap. If NULL or not provided, a default color set is used for the facet strips.
#' @param font_size A numeric value specifying the font size for the heatmap.
#' @param show_row_names A logical value indicating whether to show row names in the heatmap.
#' @param show_legend A logical value indicating whether to show the legend in the heatmap.
#' @param custom_theme A custom theme for the heatmap.
#' @param low_color A character string specifying the color for low values in the heatmap gradient. Default is "#0571b0" (blue).
#' @param mid_color A character string specifying the color for middle values in the heatmap gradient. Default is "white".
#' @param high_color A character string specifying the color for high values in the heatmap gradient. Default is "#ca0020" (red).
#' @param cluster_rows A logical value indicating whether to cluster rows (pathways). Default is FALSE.
#' @param cluster_cols A logical value indicating whether to cluster columns (samples). Default is FALSE.
#' @param clustering_method A character string specifying the clustering method. Options: "complete", "average", "single", "ward.D", "ward.D2", "mcquitty", "median", "centroid". Default is "complete".
#' @param clustering_distance A character string specifying the distance metric. Options: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "correlation", "spearman". Default is "euclidean".
#' @param dendro_line_size A numeric value specifying the line width of dendrogram branches. Default is 0.5.
#' @param dendro_labels A logical value indicating whether to show dendrogram labels. Default is FALSE.
#' @param facet_by \strong{[Deprecated]} A character string specifying an additional grouping variable for creating faceted heatmaps.
#'   This parameter is deprecated and will be removed in future versions. Use \code{secondary_groups} instead.
#' @param colorbar_title A character string specifying the title for the color bar. Default is "Z Score".
#' @param colorbar_position A character string specifying the position of the color bar. Options: "right", "left", "top", "bottom". Default is "right".
#' @param colorbar_width A numeric value specifying the width of the color bar. Default is 0.6.
#' @param colorbar_height A numeric value specifying the height of the color bar. Default is 9.
#' @param colorbar_breaks An optional numeric vector specifying custom breaks for the color bar.
#'
#' @return A ggplot heatmap object representing the heatmap of the predicted functional
#'   pathway abundance data.
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom ggh4x facet_nested strip_nested elem_list_rect
#' @importFrom stats cor as.dist dist hclust
#'
#' @examples
#' \donttest{
#' library(ggpicrust2)
#' library(ggh4x)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(magrittr)
#'
#' # Create example functional pathway abundance data
#' kegg_abundance_example <- matrix(rnorm(30), nrow = 3, ncol = 10)
#' colnames(kegg_abundance_example) <- paste0("Sample", 1:10)
#' rownames(kegg_abundance_example) <- c("PathwayA", "PathwayB", "PathwayC")
#'
#' # Create example metadata
#' metadata_example <- data.frame(
#'   sample_name = colnames(kegg_abundance_example),
#'   group = factor(rep(c("Control", "Treatment"), each = 5)),
#'   batch = factor(rep(c("Batch1", "Batch2"), times = 5))
#' )
#'
#' # Custom colors for facet strips
#' custom_colors <- c("skyblue", "salmon")
#'
#' # Example 1: Basic heatmap
#' pathway_heatmap(kegg_abundance_example, metadata_example, "group", colors = custom_colors)
#'
#' # Example 2: Heatmap with row clustering
#' pathway_heatmap(
#'   abundance = kegg_abundance_example,
#'   metadata = metadata_example,
#'   group = "group",
#'   cluster_rows = TRUE,
#'   clustering_method = "complete",
#'   clustering_distance = "euclidean",
#'   dendro_line_size = 0.8
#' )
#'
#' # Example 3: Heatmap with column clustering using correlation distance
#' pathway_heatmap(
#'   abundance = kegg_abundance_example,
#'   metadata = metadata_example,
#'   group = "group",
#'   cluster_cols = TRUE,
#'   clustering_method = "ward.D2",
#'   clustering_distance = "correlation"
#' )
#'
#' # Example 4: Multi-level grouping with secondary_groups (NEW FEATURE)
#' pathway_heatmap(
#'   abundance = kegg_abundance_example,
#'   metadata = metadata_example,
#'   group = "group",
#'   secondary_groups = "batch",
#'   colors = c("lightblue", "lightcoral", "lightgreen", "lightyellow")
#' )
#'
#' # Example 5: Custom colorbar settings
#' pathway_heatmap(
#'   abundance = kegg_abundance_example,
#'   metadata = metadata_example,
#'   group = "group",
#'   colorbar_title = "Expression Level",
#'   colorbar_position = "bottom",
#'   colorbar_width = 8,
#'   colorbar_height = 0.8,
#'   colorbar_breaks = c(-2, -1, 0, 1, 2)
#' )
#'
#' # Example 6: Advanced heatmap with clustering and custom aesthetics
#' pathway_heatmap(
#'   abundance = kegg_abundance_example,
#'   metadata = metadata_example,
#'   group = "group",
#'   cluster_rows = TRUE,
#'   cluster_cols = FALSE,  # Don't cluster columns to preserve group order
#'   clustering_method = "average",
#'   clustering_distance = "manhattan",
#'   dendro_line_size = 1.0,
#'   low_color = "#053061",     # Dark blue
#'   mid_color = "#f7f7f7",     # Light gray
#'   high_color = "#67001f",    # Dark red
#'   colorbar_title = "Z-Score",
#'   colorbar_position = "left"
#' )
#'
#' # Use real dataset
#' data("metacyc_abundance")
#' data("metadata")
#' metacyc_daa_results_df <- pathway_daa(
#'   abundance = metacyc_abundance %>% column_to_rownames("pathway"),
#'   metadata = metadata,
#'   group = "Environment",
#'   daa_method = "LinDA"
#' )
#' annotated_metacyc_daa_results_df <- pathway_annotation(
#'   pathway = "MetaCyc",
#'   daa_results_df = metacyc_daa_results_df,
#'   ko_to_kegg = FALSE
#' )
#' feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
#' 
#' # Example 7: Real data with hierarchical clustering
#' pathway_heatmap(
#'   abundance = metacyc_abundance %>%
#'     right_join(
#'       annotated_metacyc_daa_results_df %>%
#'       select(all_of(c("feature","description"))),
#'       by = c("pathway" = "feature")
#'     ) %>%
#'     filter(pathway %in% feature_with_p_0.05$feature) %>%
#'     select(-"pathway") %>%
#'     column_to_rownames("description"),
#'   metadata = metadata,
#'   group = "Environment",
#'   cluster_rows = TRUE,
#'   clustering_method = "ward.D2",
#'   clustering_distance = "correlation",
#'   colors = custom_colors,
#'   low_color = "#2166ac",  # Custom blue for low values
#'   mid_color = "#f7f7f7",  # Light gray for mid values
#'   high_color = "#b2182b", # Custom red for high values
#'   colorbar_title = "Standardized Abundance"
#' )
#'
#' # Example 8: Multiple grouping variables (NEW FEATURE)
#' # Create extended metadata with additional grouping variables
#' metadata_extended <- metadata_example %>%
#'   mutate(
#'     sex = factor(rep(c("Male", "Female"), times = 5)),
#'     age_group = factor(rep(c("Young", "Old"), each = 5))
#'   )
#'
#' # Multi-level grouping with three variables
#' pathway_heatmap(
#'   abundance = kegg_abundance_example,
#'   metadata = metadata_extended,
#'   group = "group",                    # Primary grouping
#'   secondary_groups = c("batch", "sex"), # Secondary groupings
#'   colors = c("lightblue", "lightcoral")
#' )
#'
#' # Example 9: Migration from facet_by to secondary_groups
#' # OLD WAY (deprecated, will show warning):
#' # pathway_heatmap(abundance, metadata, group = "Environment", facet_by = "Group")
#'
#' # NEW WAY (recommended):
#' # pathway_heatmap(abundance, metadata, group = "Environment", secondary_groups = "Group")
#'
#' # Example 10: Real data with multiple grouping variables
#' pathway_heatmap(
#'   abundance = metacyc_abundance %>%
#'     right_join(
#'       annotated_metacyc_daa_results_df %>%
#'       select(all_of(c("feature","description"))),
#'       by = c("pathway" = "feature")
#'     ) %>%
#'     filter(pathway %in% feature_with_p_0.05$feature) %>%
#'     select(-"pathway") %>%
#'     column_to_rownames("description"),
#'   metadata = metadata,
#'   group = "Environment",              # Primary: Pro-survival vs others
#'   secondary_groups = "Group",         # Secondary: Broad Institute vs Jackson Labs
#'   cluster_rows = TRUE,
#'   clustering_method = "ward.D2",
#'   clustering_distance = "correlation"
#' )
#' }
utils::globalVariables(c("rowname","Sample","Value","quantile","facet_nested","strip_nested","elem_list_rect", "x", "y", "xend", "yend"))

#' Generate colors for nested grouping variables
#'
#' @param metadata A data frame containing metadata
#' @param all_groups A character vector of grouping variables
#' @param colors A character vector of colors or NULL
#' @return A character vector of colors appropriate for the grouping structure
#' @keywords internal
generate_nested_colors <- function(metadata, all_groups, colors = NULL) {
  # Calculate the number of colors needed
  if (length(all_groups) == 1) {
    # Single-level grouping: one color per group level
    n_colors_needed <- length(unique(metadata[[all_groups[1]]]))
  } else {
    # Multi-level grouping: each grouping level needs distinct colors
    # Sum up unique levels across all grouping variables
    # e.g., Diet(2) + Time(4) = 6 colors needed for all facet strips
    n_colors_needed <- sum(sapply(all_groups, function(g) length(unique(metadata[[g]]))))
  }

  # Default color palette
  default_colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")

  if (is.null(colors)) {
    # Use default colors, repeat if necessary
    if (n_colors_needed <= length(default_colors)) {
      colors <- default_colors[1:n_colors_needed]
    } else {
      # Generate additional colors using colorRampPalette
      color_func <- grDevices::colorRampPalette(default_colors)
      colors <- color_func(n_colors_needed)
    }
  } else {
    # Use provided colors, repeat if necessary
    if (length(colors) < n_colors_needed) {
      warning(paste("Not enough colors provided. Need", n_colors_needed, "colors but only", length(colors), "provided. Repeating colors."))
      colors <- rep(colors, length.out = n_colors_needed)
    } else if (length(colors) > n_colors_needed) {
      colors <- colors[1:n_colors_needed]
    }
  }

  return(colors)
}

#' Create dendrogram plot from hierarchical clustering
#'
#' @param hclust_obj An hclust object from hierarchical clustering
#' @param dendro_line_size Line width for dendrogram branches
#' @param dendro_labels Whether to show labels on dendrogram
#' @param horizontal Whether to create horizontal dendrogram
#' @return A ggplot dendrogram
#' @keywords internal
create_dendrogram <- function(hclust_obj, dendro_line_size = 0.5, dendro_labels = FALSE, horizontal = FALSE) {
  if (!requireNamespace("ggdendro", quietly = TRUE)) {
    warning("Package 'ggdendro' is required for dendrogram visualization. Skipping dendrogram.")
    return(NULL)
  }
  
  # Convert hclust to dendrogram
  dendro_data <- ggdendro::dendro_data(hclust_obj)
  
  # Create the plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = dendro_data$segments, 
                         ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                         linewidth = dendro_line_size, color = "black") +
    ggplot2::theme_void()
  
  if (dendro_labels) {
    p <- p + ggplot2::geom_text(data = dendro_data$labels, 
                               ggplot2::aes(x = x, y = y, label = label),
                               size = 3, hjust = 0.5, vjust = 1)
  }
  
  if (horizontal) {
    p <- p + ggplot2::coord_flip()
  }
  
  return(p)
}
pathway_heatmap <- function(abundance,
                            metadata,
                            group,
                            secondary_groups = NULL,
                            colors = NULL,
                            font_size = 12,
                            show_row_names = TRUE,
                            show_legend = TRUE,
                            custom_theme = NULL,
                            low_color = "#0571b0",
                            mid_color = "white",
                            high_color = "#ca0020",
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            clustering_method = "complete",
                            clustering_distance = "euclidean",
                            dendro_line_size = 0.5,
                            dendro_labels = FALSE,
                            facet_by = NULL,
                            colorbar_title = "Z Score",
                            colorbar_position = "right",
                            colorbar_width = 0.6,
                            colorbar_height = 9,
                            colorbar_breaks = NULL) {
  # Input validation using unified functions
  validate_abundance(abundance, min_samples = 2)
  validate_metadata(metadata)

  # Ensure abundance is a matrix with names
  abundance <- as.matrix(abundance)
  if (nrow(abundance) < 1) stop("At least one pathway is required")
  if (is.null(colnames(abundance))) colnames(abundance) <- paste0("Sample", seq_len(ncol(abundance)))
  if (is.null(rownames(abundance))) rownames(abundance) <- paste0("Pathway", seq_len(nrow(abundance)))

  # Handle deprecated facet_by parameter
 if (!is.null(facet_by)) {
    warning("'facet_by' is deprecated. Use 'secondary_groups' instead.", call. = FALSE)
    if (is.null(secondary_groups)) secondary_groups <- facet_by
  }

  # Build and validate all grouping variables
  all_groups <- c(group, secondary_groups)
  for (grp in all_groups) {
    validate_group(metadata, grp, min_groups = 2)
  }

  if (!is.null(colors) && !is.character(colors)) {
    stop("colors must be NULL or a character vector")
  }
  
  # Heatmaps use color changes to visualize changes in values. However, if the
  # data for plotting the heat map are too different, for example, if the heat
  # map is plotted using gene expression data, gene1 is expressed above 1000 in
  # all samples and gene2 is expressed between 1-10 in all samples, it is
  # difficult to plot the heat map with small changes in the expression of two
  # genes in different samples by the colors to reflect. Therefore, when
  # plotting a heat map, we usually normalize the gene expression data, that
  # is, we subtract the mean value of each gene expression from the expression
  # of this gene in all samples and divide it by its standard deviation, and
  # this normalization is called standard normalization or Z-score processing.
  # The processed values are reduced equally, and the expression of each gene
  # in all samples becomes a set of values with a mean of 0 and a standard
  # deviation of 1. At this point, the plotted heat map gives a good indication
  # of the variation in expression of all genes across samples.

  # Align samples using unified function
  aligned <- align_samples(abundance, metadata)
  abundance <- as.matrix(aligned$abundance)
  metadata <- aligned$metadata
  metadata$sample_name <- colnames(abundance)

  # Validate group column
  validate_group(metadata, group, min_groups = 2)

  # Perform z-score normalization
  z_abundance <- t(apply(abundance, 1, scale))
  colnames(z_abundance) <- colnames(abundance)

  # Convert the abundance matrix to a data frame
  z_df <- as.data.frame(z_abundance)

  metadata <- metadata %>% as.data.frame()

  # Perform clustering if requested
  row_order <- rownames(z_abundance)
  col_order <- colnames(z_abundance)
  
  if (cluster_rows) {
    # Calculate distance matrix for rows
    if (clustering_distance == "correlation") {
      # Handle correlation distance specially
      cor_matrix <- cor(t(z_abundance))
      row_dist <- as.dist(1 - cor_matrix)
    } else if (clustering_distance == "spearman") {
      # Handle spearman correlation distance specially
      cor_matrix <- cor(t(z_abundance), method = "spearman")
      row_dist <- as.dist(1 - cor_matrix)
    } else {
      row_dist <- dist(z_abundance, method = clustering_distance)
    }
    row_hclust <- hclust(row_dist, method = clustering_method)
    row_order <- rownames(z_abundance)[row_hclust$order]
  }
  
  # Always prepare ordered metadata for consistent behavior
  # Order by all grouping variables hierarchically
  if (length(all_groups) == 1) {
    ordered_metadata <- metadata[order(metadata[, all_groups[1]]),]
  } else {
    # Order by multiple grouping variables
    ordered_metadata <- metadata[do.call(order, metadata[all_groups]),]
  }
  ordered_sample_names <- ordered_metadata$sample_name
  ordered_group_levels <- ordered_metadata %>% select(all_of(c(group))) %>% pull()
  
  if (cluster_cols) {
    # Calculate distance matrix for columns
    if (clustering_distance == "correlation") {
      # Handle correlation distance specially
      cor_matrix <- cor(z_abundance)
      col_dist <- as.dist(1 - cor_matrix)
    } else if (clustering_distance == "spearman") {
      # Handle spearman correlation distance specially
      cor_matrix <- cor(z_abundance, method = "spearman")
      col_dist <- as.dist(1 - cor_matrix)
    } else {
      col_dist <- dist(t(z_abundance), method = clustering_distance)
    }
    col_hclust <- hclust(col_dist, method = clustering_method)
    col_order <- colnames(z_abundance)[col_hclust$order]
  } else {
    # Order the samples based on the environment information (default behavior)
    col_order <- ordered_sample_names
  }

  # Convert the abundance data frame to a long format
  # Prepare metadata columns to join - include all grouping variables
  metadata_cols <- c("sample_name", all_groups)
  
  long_df <- z_df %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(cols = -rowname,
                        names_to = "Sample",
                        values_to = "Value") %>% 
    left_join(metadata %>% select(all_of(metadata_cols)), 
              by = c("Sample" = "sample_name"))

  # Set the order of the samples and pathways in the heatmap
  long_df$Sample <- factor(long_df$Sample, levels = col_order)
  long_df$rowname <- factor(long_df$rowname, levels = row_order)

  # Compute breaks from the data
  breaks <- range(long_df$Value, na.rm = TRUE)

  # Generate appropriate colors for the grouping structure
  colors <- generate_nested_colors(metadata, all_groups, colors)

  # Create the heatmap using ggplot
  p <- ggplot2::ggplot(data = long_df,
                    mapping = ggplot2::aes(x = Sample, y = rowname, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = 0) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_y_discrete(expand = c(0, 0), position = "left") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    # Customize the appearance of the heatmap
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = font_size, color = "black"),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        color = "black",
        size = 10,
        face = "bold"
      ),
      panel.spacing = unit(0, "lines"),
      legend.title = ggplot2::element_text(size = 12, color = "black",face = "bold"),
      legend.text = ggplot2::element_text(size = 12, color = "black",face = "bold"),
      panel.background = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(l = 0, unit = "cm"),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    # Add improved color bar to the heatmap
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        direction = if(colorbar_position %in% c("left", "right")) "vertical" else "horizontal",
        reverse = FALSE,
        barwidth = unit(colorbar_width, "cm"),
        barheight = unit(colorbar_height, "cm"),
        title = colorbar_title,
        title.position = if(colorbar_position %in% c("left", "right")) "top" else "left",
        title.hjust = if(colorbar_position %in% c("left", "right")) 0.5 else 0,
        ticks = TRUE,
        label = TRUE,
        breaks = colorbar_breaks
      )
    ) + 
    ggplot2::theme(legend.position = colorbar_position)

  # Add faceting support with multiple grouping levels
  if (length(all_groups) == 1) {
    # Single-level faceting
    p <- p + ggh4x::facet_nested(
      cols = vars(!!sym(all_groups[1])),
      space = "free",
      scale = "free",
      switch = "x",
      strip = ggh4x::strip_nested(
        background_x = ggh4x::elem_list_rect(fill = colors)
      )
    )
  } else {
    # Multi-level nested faceting
    p <- p + ggh4x::facet_nested(
      cols = vars(!!!syms(all_groups)),
      space = "free",
      scale = "free",
      switch = "x",
      strip = ggh4x::strip_nested(
        background_x = ggh4x::elem_list_rect(fill = colors)
      )
    )
  }

  if (!show_row_names) {
    p <- p + theme(axis.text.y = element_blank())
  }
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (!is.null(custom_theme)) {
    p <- p + custom_theme
  }

  # Create dendrograms if clustering was performed
  if (cluster_rows || cluster_cols) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      warning("Package 'patchwork' is required for combining plots with dendrograms. Returning heatmap only.")
      return(p)
    }
    
    # Create row dendrogram if rows were clustered
    if (cluster_rows && exists("row_hclust")) {
      row_dendro <- create_dendrogram(row_hclust, dendro_line_size, dendro_labels, horizontal = TRUE)
      if (!is.null(row_dendro)) {
        p <- row_dendro + p + patchwork::plot_layout(widths = c(0.2, 1))
      }
    }
    
    # Create column dendrogram if columns were clustered
    if (cluster_cols && exists("col_hclust")) {
      col_dendro <- create_dendrogram(col_hclust, dendro_line_size, dendro_labels, horizontal = FALSE)
      if (!is.null(col_dendro)) {
        p <- col_dendro / p + patchwork::plot_layout(heights = c(0.2, 1))
      }
    }
  }

  # Print the ordered sample names and group levels
  if (!cluster_cols) {
    cat("The Sample Names in order from left to right are:\n")
    cat(ordered_sample_names, sep = ", ")
    cat("\n")

    cat("The Group Levels in order from left to right are:\n")
    cat(ordered_group_levels, sep = ", ")
    cat("\n")
  } else {
    cat("Samples ordered by hierarchical clustering (", clustering_method, " method, ", clustering_distance, " distance)\n")
    cat("Column order: ", paste(col_order, collapse = ", "), "\n")
  }
  
  if (cluster_rows) {
    cat("Pathways ordered by hierarchical clustering (", clustering_method, " method, ", clustering_distance, " distance)\n")
    cat("Row order: ", paste(row_order, collapse = ", "), "\n")
  }

  return(p)
}
