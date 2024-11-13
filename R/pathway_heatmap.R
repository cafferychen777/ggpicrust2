#' Create pathway heatmap
#'
#' This function creates a heatmap of the predicted functional pathway abundance data.
#' The function first makes the abundance data relative, then converts the abundance data
#' to a long format and orders the samples based on the environment information.
#' The heatmap is then created using the `ggplot2` library.
#'
#' @name pathway_heatmap
#' @param abundance A matrix or data frame of pathway abundance data, where columns
#'   correspond to samples and rows correspond to pathways. Must contain at least
#'   two samples.
#' @param metadata A data frame of metadata, where each row corresponds to a sample
#'   and each column corresponds to a metadata variable.
#' @param group A character string specifying the column name in the metadata data frame
#'   that contains the group variable. Must contain at least two groups.
#' @param colors A vector of colors used for the background of the facet labels in the
#'   heatmap. If NULL or not provided, a default color set is used for the facet strips.
#' @param font_size A numeric value specifying the font size for the heatmap.
#' @param show_row_names A logical value indicating whether to show row names in the heatmap.
#' @param show_legend A logical value indicating whether to show the legend in the heatmap.
#' @param custom_theme A custom theme for the heatmap.
#'
#' @return A ggplot heatmap object representing the heatmap of the predicted functional
#'   pathway abundance data.
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom ggh4x facet_nested strip_nested elem_list_rect
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
#'   group = factor(rep(c("Control", "Treatment"), each = 5))
#' )
#'
#' # Custom colors for facet strips
#' custom_colors <- c("skyblue", "salmon")
#'
#' # Create a heatmap using custom colors for facet strips
#' pathway_heatmap(kegg_abundance_example, metadata_example, "group", colors = custom_colors)
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
#'   colors = custom_colors
#' )
#' }
utils::globalVariables(c("rowname","Sample","Value","quantile","facet_nested","strip_nested","elem_list_rect"))
pathway_heatmap <- function(abundance,
                            metadata,
                            group,
                            colors = NULL,
                            font_size = 12,
                            show_row_names = TRUE,
                            show_legend = TRUE,
                            custom_theme = NULL) {
  # Input validation
  if (!is.matrix(abundance) && !is.data.frame(abundance)) {
    stop("abundance must be a data frame or matrix")
  }
  
  # Ensure abundance is a matrix
  abundance <- as.matrix(abundance)
  
  # Check sample count
  if (ncol(abundance) < 2) {
    stop("At least two samples are required for creating a heatmap")
  }
  
  # Check pathway count
  if (nrow(abundance) < 1) {
    stop("At least one pathway is required")
  }
  
  # Ensure column names exist
  if (is.null(colnames(abundance))) {
    colnames(abundance) <- paste0("Sample", seq_len(ncol(abundance)))
  }
  
  # Ensure row names exist
  if (is.null(rownames(abundance))) {
    rownames(abundance) <- paste0("Pathway", seq_len(nrow(abundance)))
  }
  
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  
  if (!is.character(group) || length(group) != 1) {
    stop("group must be a single character string")
  }
  
  if (!is.null(colors) && !is.character(colors)) {
    stop("colors must be NULL or a character vector of color codes")
  }
  
  # Check group count
  group_levels <- unique(metadata[[group]])
  if (length(group_levels) < 2) {
    stop("At least two groups are required for comparison")
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

  # Check that 'group' is a column in 'metadata'
  if (!group %in% colnames(metadata)) {
    stop(paste("group:", group, "must be a column in metadata"))
  }

  # Find the column in metadata that matches the column names of abundance
  sample_name_col <- colnames(metadata)[sapply(colnames(metadata), function(x) all(colnames(abundance) %in% metadata[[x]]))]
  metadata$sample_name <- metadata %>% select(all_of(c(sample_name_col))) %>% pull()

  if (!all(colnames(abundance) %in% metadata$sample_name)) {
    stop("Samples in abundance and metadata must match")
  }

  # Perform z-score normalization
  z_abundance <- t(apply(abundance, 1, scale))
  colnames(z_abundance) <- colnames(abundance)

  # Convert the abundance matrix to a data frame
  z_df <- as.data.frame(z_abundance)

  metadata <- metadata %>% as.data.frame()

  # Order the samples based on the environment information
  ordered_metadata <- metadata[order(metadata[, group]),]
  ordered_sample_names <- ordered_metadata$sample_name
  order <- ordered_metadata$sample_name
  ordered_group_levels <- ordered_metadata %>% select(all_of(c(group))) %>% pull()

  # Convert the abundance data frame to a long format
  long_df <- z_df %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(cols = -rowname,
                        names_to = "Sample",
                        values_to = "Value") %>% 
    left_join(metadata %>% select(all_of(c("sample_name",group))), 
              by = c("Sample" = "sample_name"))

  # Set the order of the samples in the heatmap
  long_df$Sample <- factor(long_df$Sample, levels = order)

  # Compute breaks from the data
  breaks <- range(long_df$Value, na.rm = TRUE)

  if (is.null(colors)) {
    colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")
  }

  # Create the heatmap using ggplot
  p <- ggplot2::ggplot(data = long_df,
                    mapping = ggplot2::aes(x = Sample, y = rowname, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020", midpoint = 0) +
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
    # Add a color bar to the heatmap
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        direction = "vertical",
        reverse = F,
        barwidth = unit(0.6, "cm"),
        barheight = unit(9, "cm"),
        title = "Z Score",
        title.position = "top",
        title.hjust = -1,
        ticks = TRUE,
        label = TRUE
      )
    ) + 
    ggh4x::facet_nested(
      cols = vars(!!sym(group)), 
      space = "free", 
      scale = "free", 
      switch = "x", 
      strip = ggh4x::strip_nested(
        background_x = ggh4x::elem_list_rect(fill = colors)
      )
    )

  if (!show_row_names) {
    p <- p + theme(axis.text.y = element_blank())
  }
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (!is.null(custom_theme)) {
    p <- p + custom_theme
  }

  # Print the ordered sample names and group levels
  cat("The Sample Names in order from left to right are:\n")
  cat(ordered_sample_names, sep = ", ")
  cat("\n")

  cat("The Group Levels in order from left to right are:\n")
  cat(ordered_group_levels, sep = ", ")
  cat("\n")

  return(p)
}
